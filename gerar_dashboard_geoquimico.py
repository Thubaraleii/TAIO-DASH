"""Dashboard geoquimico -- terceiro produto web (lista + mapa + diagramas),
100% escopo Taio (a versao anterior misturava com uma planilha de geoquimica
da PMP/Florianopolis -- removida por pedido explicito do usuario; o produto
PMP vai virar um dashboard SEPARADO, sem afetar este aqui).

Dados: 308 pontos de campo de `pontos_unificados_completo.gpkg` (mesmo
catalogo unificado usado nos outros 3 produtos), sem geoquimica bruta -- so
um flag Sim/Nao + classificacao Alto/Baixo-Ti pra ~23 pontos.

Os diagramas geoquimicos (TAS, AFM, Shand's Index, e 3 campos de magma-tipo
digitalizados de figura da literatura) mostram SO os campos de referencia,
sem pontos plotados -- o Taio ainda nao tem oxido/traço bruto de verdade pra
posicionar amostras neles. Quando existir esse dado, e so adicionar os
pontos nas mesmas figuras.

Mapa: Leaflet (tiles reais, zoom/pan continuo) -- mesma tecnica do webmap
dedicado (gerar_webmap_taio.py), NAO mais Plotly com imagem estatica
esticada (desempenho/nitidez bem inferior, e tinha um bug real de
renderizacao: Plotly deduplica layout images por posicao, nao indice).

Uso:
    python visualizacao_web/gerar_dashboard_geoquimico.py

Gera:
    visualizacao_web/dashboard_geoquimico.html
"""
import base64
import json
from collections import Counter
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from pyproj import Transformer

BASE = Path(__file__).parent
PONTOS_CAMPO_GPKG = (
    BASE.parent.parent / "2_Banco_de_Dados" / "Unificação" / "GPKG_Novos" / "pontos_unificados_completo.gpkg"
)
LOGO_PATH = BASE / "assets" / "logo_gstech.jpg"
OUT_HTML = BASE / "dashboard_geoquimico.html"

# identidade visual GS Tech -- mesma paleta dos outros produtos, manter em sincronia.
MARCA_ROXO_ESCURO = "#2D0A4A"
MARCA_ROXO = "#7B2FFF"
MARCA_AZUL = "#2E6F95"
MARCA_NAVY = "#1B1F2E"
MARCA_CINZA_CLARO = "#F2F2F2"
MARCA_FONTE = "Montserrat, Arial, sans-serif"
COR_PAINEL = "#262B3D"

COR_SILL = "#A63D2F"
COR_DIQUE = "#1B4332"
CORES_LITOLOGIA_CAMPO = {
    "sill_diabasio": COR_SILL, "sill_diabasio_cprm": COR_SILL,
    "dique": COR_DIQUE, "dique_cprm": COR_DIQUE,
    "encaixante_teresina": "#D6C79A", "encaixante_serra_alta": "#8C8C86",
    "encaixante_irati": "#3E362C", "encaixante_palermo": "#B5AE93",
    "encaixante_rio_bonito": "#C9A66B", "encaixante_sedimentar": "#999999",
}
COR_LITOLOGIA_PADRAO = "#999999"
COR_TI_INDEFINIDO = "#9A9A9A"


def txt(valor, padrao=""):
    """Converte valor de celula (pode ser pd.NA/None/nan) pra string segura --
    `pd.NA or ""` explode com TypeError, entao nao da pra usar `or` direto."""
    if valor is None or (isinstance(valor, float) and np.isnan(valor)) or valor is pd.NA:
        return padrao
    try:
        if pd.isna(valor):
            return padrao
    except (TypeError, ValueError):
        pass
    return str(valor)


def logo_base64():
    if not LOGO_PATH.exists():
        return None
    return base64.b64encode(LOGO_PATH.read_bytes()).decode("ascii")


# ======================================================================
# 1. dados de campo (Taio) -- unica fonte deste dashboard
# ======================================================================
def carregar_campo():
    gdf = gpd.read_file(PONTOS_CAMPO_GPKG)
    if gdf.crs is None or gdf.crs.to_epsg() != 31982:
        gdf = gdf.to_crs(31982)
    registros = []
    for row in gdf.itertuples():
        desc = row.descricao_campo if pd.notna(row.descricao_campo) else ""
        registros.append({
            "id": f"CAMPO-{row.ponto_id}",
            "nome": row.ponto_id,
            "litologia": row.litologia_padronizada if pd.notna(row.litologia_padronizada) else "indefinido",
            "classificacao_ti": row.ti_geoquimico if pd.notna(row.ti_geoquimico) else None,
            "tipo_ponto": row.tipo_ponto if pd.notna(row.tipo_ponto) else "",
            "qualidade": row.qualidade_dado if pd.notna(row.qualidade_dado) else "",
            "x": float(row.X) if pd.notna(row.X) else None,
            "y": float(row.Y) if pd.notna(row.Y) else None,
            "descricao": desc,
            "tem_geoquimica": row.geoquimica == "Sim",
            "cor_mapa": CORES_LITOLOGIA_CAMPO.get(row.litologia_padronizada, COR_LITOLOGIA_PADRAO),
        })
    return registros


# ======================================================================
# 2. suavizacao de poligonos (campos de literatura) -- Catmull-Rom, mesma
#    tecnica portada do notebook GeoQMASTER_VF.ipynb do usuario.
# ======================================================================
def suavizar_poligono_catmull_rom(pontos, n_interpolacao=25):
    pontos = np.array(pontos, dtype=float)
    pts = np.vstack([pontos[-2], pontos[-1], pontos, pontos[0], pontos[1]])
    curva = []
    for i in range(2, len(pts) - 2):
        p0, p1, p2, p3 = pts[i - 1], pts[i], pts[i + 1], pts[i + 2]
        for t in np.linspace(0, 1, n_interpolacao):
            t2, t3 = t * t, t * t * t
            ponto = 0.5 * (
                (2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2
                + (-p0 + 3 * p1 - 3 * p2 + p3) * t3
            )
            curva.append(ponto)
    return np.array(curva)


def montar_shapes_campos(campos, n_interpolacao=25, suavizar=True):
    """Devolve (lista de go.Scatter preenchidos, lista de anotacoes) pros
    campos de literatura de um diagrama. suavizar=False desenha o poligono
    cru (usado pros retangulos da Tabela 1 -- arestas retas comunicam
    "intervalo numerico publicado", nao um contorno organico digitalizado)."""
    traces, anotacoes = [], []
    for nome, info in campos.items():
        if suavizar:
            curva = suavizar_poligono_catmull_rom(info["coords"], n_interpolacao=n_interpolacao)
            xs, ys = curva[:, 0].tolist(), curva[:, 1].tolist()
        else:
            xs, ys = [p[0] for p in info["coords"]], [p[1] for p in info["coords"]]
        traces.append(go.Scatter(
            x=xs + [xs[0]], y=ys + [ys[0]],
            mode="lines", fill="toself", fillcolor=info["cor"], opacity=0.55 if not suavizar else 0.72,
            line=dict(color="black", width=1.2), hoverinfo="skip", showlegend=False,
        ))
        anotacoes.append(dict(
            x=info["label"][0], y=info["label"][1], xref="x", yref="y",
            text=f"<b>{nome}</b>", showarrow=False, font=dict(size=10 if not suavizar else 11, color="black"),
        ))
    return traces, anotacoes


def tema_grafico(fig, titulo, altura=340, nota=None):
    fig.update_layout(
        paper_bgcolor=MARCA_NAVY, plot_bgcolor=COR_PAINEL, height=altura + (26 if nota else 0),
        font=dict(family=MARCA_FONTE, color=MARCA_CINZA_CLARO, size=11),
        margin=dict(l=55, r=20, t=45, b=45 + (26 if nota else 0)),
        title=dict(text=f"<b>{titulo}</b>", x=0.02, font=dict(size=13, color=MARCA_CINZA_CLARO)),
        legend=dict(bgcolor="rgba(45,10,74,0.75)", bordercolor=MARCA_ROXO, borderwidth=1, font=dict(size=10)),
    )
    fig.update_xaxes(color=MARCA_CINZA_CLARO, gridcolor="#3a3f52", zerolinecolor="#3a3f52")
    fig.update_yaxes(color=MARCA_CINZA_CLARO, gridcolor="#3a3f52", zerolinecolor="#3a3f52")
    if nota:
        # posicionado em coordenadas de PIXEL (nao "paper") ancorado no canto
        # inferior direito, abaixo do titulo do eixo X -- em "paper" a nota
        # colidia com o titulo do eixo (ambos ficavam na mesma faixa vertical
        # negativa, sobrepondo o texto).
        fig.add_annotation(
            text=nota, xref="paper", yref="paper", x=1, y=0, xshift=0, yshift=-58,
            showarrow=False, font=dict(size=9, color=MARCA_CINZA_CLARO), opacity=0.55, xanchor="right", yanchor="top",
        )
    return fig


NOTA_SEM_DADO = "Só campos de referência — Taió ainda não tem óxido/traço bruto pra plotar"
NOTA_APROXIMADO = NOTA_SEM_DADO + " · campos digitalizados aproximados (Fontoura, TCC 2024, Fig.12/33 · Peate et al. 1997)"


# ======================================================================
# 3. diagramas geoquimicos padrao (TAS, AFM, Shand's Index) -- limites bem
#    conhecidos da literatura (Le Bas et al. 1986 / Irvine & Baragar 1971),
#    simplificados (linhas principais, sem todas as subdivisoes finas) ja
#    que ainda nao ha ponto real do Taio pra classificar contra eles.
# ======================================================================
def montar_tas():
    fig = go.Figure()
    linha = dict(color="#A63D2F", width=1.3)
    # curva alcalino/sub-alcalino (aproximada, tracejada)
    curva_x = [39, 41, 45, 50, 55, 60, 65, 70, 75, 78]
    curva_y = [0.3, 0.6, 1.2, 2.5, 4.0, 5.7, 7.3, 8.6, 9.6, 10.2]
    fig.add_trace(go.Scatter(x=curva_x, y=curva_y, mode="lines", line=dict(color="#A63D2F", width=1.2, dash="dot"),
                              hoverinfo="skip", showlegend=False))
    # divisorias verticais principais (base da serie sub-alcalina)
    for sio2, ymax in [(45, 3), (52, 5), (57, 5.9), (63, 7), (69, 8)]:
        fig.add_trace(go.Scatter(x=[sio2, sio2], y=[0, ymax], mode="lines", line=linha,
                                  hoverinfo="skip", showlegend=False))
    # contorno externo aproximado (moldura simplificada dos campos alcalinos)
    contorno = [
        (41, 3), (41, 7), (45, 9.5), (48.5, 11.5), (52.5, 14), (57, 14.2), (61, 11.8),
        (63, 9), (69, 8.6), (69.5, 12.5), (73, 14.5), (77, 10.2), (78, 6), (69, 3.2),
        (63, 0), (45, 0), (41, 3),
    ]
    fig.add_trace(go.Scatter(x=[p[0] for p in contorno], y=[p[1] for p in contorno], mode="lines",
                              line=linha, hoverinfo="skip", showlegend=False))
    rotulos = [
        ("Picrobasalto", 43, 0.8), ("Basalto", 48.5, 1.7), ("Andesito\nbasáltico", 54.5, 2.2),
        ("Andesito", 60, 2.5), ("Dacito", 66, 2.8), ("Riolito", 74, 5),
        ("Traqui-\nbasalto", 47.5, 4.6), ("Traquiandesito\nbasáltico", 54, 5.6),
        ("Traqui-\nandesito", 59.5, 7), ("Traquito\nTraquidacito", 66, 10),
        ("Basanito\nTefrítico", 43, 5.5), ("Fono-\nTefrito", 48, 8), ("Tefrito\nFonolítico", 53.5, 11.3),
        ("Fonolito", 55, 14), ("Foiditos", 40.5, 5), ("Alcalino", 41, 1.2), ("Sub-alcalino/Toleiítico", 50, 0.4),
    ]
    for texto, x, y in rotulos:
        fig.add_annotation(x=x, y=y, xref="x", yref="y", text=texto.replace("\n", "<br>"), showarrow=False,
                            font=dict(size=8.5, color="#A63D2F"), align="center")
    fig.update_xaxes(title_text="SiO₂", range=[38, 80])
    fig.update_yaxes(title_text="Na₂O + K₂O", range=[0, 16])
    return tema_grafico(fig, "TAS — Total Álcalis vs Sílica", altura=380, nota=NOTA_SEM_DADO)


def montar_afm():
    fig = go.Figure()
    # curva de Irvine & Baragar (1971), aproximada (A=alcalis, F=FeOt, M=MgO;
    # Plotly ternario: 'a' = vertice de cima, 'b' = esquerda, 'c' = direita
    # -- aqui a=F(topo), b=A(esquerda), c=M(direita), igual a imagem de referencia)
    curva = [
        (61, 0, 39), (57, 8, 35), (50, 15, 35), (44, 22, 34), (41, 30, 29),
        (37, 42, 21), (33, 55, 12), (30, 68, 2),
    ]
    fig.add_trace(go.Scatterternary(
        a=[p[0] for p in curva], b=[p[1] for p in curva], c=[p[2] for p in curva],
        mode="lines", line=dict(color="#A63D2F", width=1.4), hoverinfo="skip", showlegend=False,
    ))
    fig.update_layout(
        ternary=dict(
            sum=100, bgcolor=COR_PAINEL,
            aaxis=dict(title="F", color=MARCA_CINZA_CLARO, gridcolor="#3a3f52", linecolor=MARCA_CINZA_CLARO),
            baxis=dict(title="A", color=MARCA_CINZA_CLARO, gridcolor="#3a3f52", linecolor=MARCA_CINZA_CLARO),
            caxis=dict(title="M", color=MARCA_CINZA_CLARO, gridcolor="#3a3f52", linecolor=MARCA_CINZA_CLARO),
        ),
        annotations=[
            dict(text="Toleítico", x=0.5, y=0.62, xref="paper", yref="paper", showarrow=False,
                 font=dict(size=11, color="#A63D2F")),
            dict(text="Cálcio-alcalino", x=0.5, y=0.18, xref="paper", yref="paper", showarrow=False,
                 font=dict(size=11, color="#A63D2F")),
        ],
    )
    return tema_grafico(fig, "AFM — Álcalis / FeOt / MgO", altura=380, nota=NOTA_SEM_DADO)


def montar_shand():
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[1, 1], y=[0, 7], mode="lines", line=dict(color="#A63D2F", width=1.3),
                              hoverinfo="skip", showlegend=False))
    fig.add_trace(go.Scatter(x=[0, 1.6], y=[1, 1], mode="lines", line=dict(color="#A63D2F", width=1.3),
                              hoverinfo="skip", showlegend=False))
    fig.add_trace(go.Scatter(x=[0.4, 1.5], y=[0.4, 1.65], mode="lines",
                              line=dict(color="#A63D2F", width=1.2, dash="dot"), hoverinfo="skip", showlegend=False))
    for texto, x, y in [
        ("Metaluminoso", 0.65, 6), ("Peraluminoso", 1.3, 6),
        ("Peralcalino", 0.55, 0.35),
    ]:
        fig.add_annotation(x=x, y=y, xref="x", yref="y", text=texto, showarrow=False,
                            font=dict(size=10, color="#A63D2F"))
    fig.update_xaxes(title_text="A/CNK", range=[0, 1.6])
    fig.update_yaxes(title_text="A/NK", range=[0, 7])
    return tema_grafico(fig, "Shand's Index — A/NK × A/CNK", nota=NOTA_SEM_DADO)


# ======================================================================
# 4. campos de magma-tipo -- fonte: Fontoura, G.M. (TCC UFSC, 2024),
#    "Geologia, petrografia e geoquímica elemental dos diques máficos do
#    Farol de Santa Marta, SC" (Figura 12, p.46, "Modificado de Peate et
#    al., 1997") e Tabela 1 (p.43, "Peate et al., 1992").
#
#    Painel I (MgO×TiO2, Ti/Y×Ti/Zr) no original mostra campos suaves
#    (blob) -- vertices digitalizados a olho na imagem extraida do PDF
#    (melhor resolucao que o print enviado, mas ainda aproximado, nao sao
#    coordenadas numericas publicadas).
#
#    Paineis II/III (Fe2O3×TiO2, Zr/Y×Sr, Zr/Y×Ti/Zr, Sr×Ti/Y) no original
#    mostram NUVENS DE PONTOS (amostras de referencia), nao campos --
#    como nao ha como digitalizar pontos individuais com precisao, uso os
#    INTERVALOS NUMERICOS PUBLICADOS na Tabela 1 (Peate et al. 1992) como
#    retangulo de referencia (min/max de cada eixo por magma-tipo) -- mais
#    fiel ao dado real citavel do que tentar desenhar a nuvem por olho.
# ======================================================================
CAMPOS_MGO_TIO2 = {
    "Urubici": {"cor": "#B4D900", "label": (3.75, 4.3), "coords": [
        (3.0, 3.3), (3.1, 3.9), (3.4, 4.3), (3.9, 4.5), (4.3, 4.3), (4.5, 3.9),
        (4.3, 3.4), (3.9, 3.15), (3.4, 3.1),
    ]},
    "Pitanga": {"cor": "#E6951A", "label": (4.35, 3.6), "coords": [
        (3.5, 2.9), (3.6, 3.4), (3.9, 3.8), (4.4, 4.0), (4.9, 3.85), (5.2, 3.4),
        (5.0, 2.95), (4.5, 2.7), (4.0, 2.65),
    ]},
    "Paranapanema": {"cor": "#B9D7F0", "label": (4.6, 2.05), "coords": [
        (2.8, 1.75), (2.9, 2.1), (3.3, 2.45), (4.0, 2.6), (4.8, 2.55), (5.5, 2.35),
        (6.1, 2.05), (6.4, 1.75), (6.0, 1.55), (5.2, 1.5), (4.2, 1.55), (3.4, 1.6), (2.9, 1.65),
    ]},
    "Ribeira": {"cor": "#8791D0", "label": (2.9, 2.1), "coords": [
        (2.55, 2.0), (2.6, 2.25), (2.9, 2.35), (3.25, 2.25), (3.3, 2.0), (3.0, 1.9), (2.7, 1.9),
    ]},
    "Esmeralda": {"cor": "#10962F", "label": (5.6, 1.4), "coords": [
        (4.2, 1.15), (4.4, 1.55), (4.9, 1.85), (5.6, 1.95), (6.3, 1.8), (6.8, 1.5),
        (6.9, 1.15), (6.5, 0.95), (5.8, 0.9), (5.0, 0.9), (4.5, 0.95),
    ]},
    "Gramado": {"cor": "#D62728", "label": (5.5, 0.95), "coords": [
        (2.6, 0.85), (2.65, 1.15), (2.9, 1.4), (3.4, 1.55), (4.0, 1.5), (4.6, 1.35),
        (5.2, 1.25), (5.8, 1.2), (6.5, 1.15), (7.3, 1.05), (8.0, 0.95), (8.6, 0.85),
        (8.8, 0.7), (8.3, 0.55), (7.3, 0.5), (6.0, 0.5), (4.8, 0.55), (3.7, 0.6), (2.9, 0.65),
    ]},
}
CAMPOS_TIY_TIZR = {
    "Ribeira": {"cor": "#8791D0", "label": (282, 87), "coords": [
        (260, 82), (270, 90), (285, 92), (300, 88), (305, 80), (295, 73), (278, 72), (265, 76),
    ]},
    "Paranapanema": {"cor": "#B9D7F0", "label": (355, 92), "coords": [
        (300, 65), (305, 80), (315, 90), (335, 96), (360, 97), (385, 93), (400, 85),
        (405, 75), (395, 68), (370, 64), (340, 62), (315, 62),
    ]},
    "Pitanga": {"cor": "#E6951A", "label": (520, 90), "coords": [
        (400, 66), (410, 80), (430, 90), (460, 96), (500, 98), (550, 96), (600, 90),
        (635, 80), (645, 70), (630, 65), (590, 63), (540, 62), (480, 62), (430, 63),
    ]},
    "Urubici": {"cor": "#B4D900", "label": (580, 70), "coords": [
        (480, 60), (490, 68), (510, 74), (540, 76), (580, 75), (620, 72), (655, 68),
        (670, 63), (655, 58), (610, 56), (560, 56), (510, 57),
    ]},
    "Esmeralda": {"cor": "#10962F", "label": (225, 70), "coords": [
        (195, 58), (200, 68), (220, 76), (245, 78), (258, 72), (255, 63), (240, 58), (215, 55),
    ]},
    "Gramado": {"cor": "#D62728", "label": (205, 45), "coords": [
        (175, 37), (180, 45), (195, 52), (215, 57), (230, 55), (235, 48), (225, 40), (205, 35), (185, 33),
    ]},
}

# Paineis C/D/E/F -- digitalizados da FIGURA 33 (p.89 do TCC, "Fonte: Do
# autor, 2024"), nao mais da Figura 12/Tabela 1. A Fig.33 e a reproducao do
# proprio autor com os MESMOS campos suaves (organicos) da Fig.12, so que
# com as amostras dele plotadas por cima -- exportei a imagem embutida no
# PDF em resolucao alta (bem melhor que o print originalmente enviado) e
# digitalizei os vertices dos campos direto dela. Ainda e aproximado (nao
# sao coordenadas numericas publicadas), mas com uma fonte de imagem bem
# mais nitida que antes.
CAMPOS_C_FE2O3_TIO2 = {
    "Paranapanema": {"cor": "#B5B5B5", "label": (2.2, 15.5), "coords": [
        (1.5, 13), (1.6, 14.5), (1.9, 16.2), (2.3, 17.5), (2.8, 17.3), (3.1, 16),
        (3.2, 14.5), (3.0, 13), (2.5, 12.3), (1.9, 12.5),
    ]},
    "Pitanga": {"cor": "#E6951A", "label": (3.5, 15), "coords": [
        (2.9, 13), (3.0, 14.5), (3.3, 16), (3.7, 17.3), (4.0, 16.8), (4.1, 15),
        (4.0, 13.3), (3.6, 12.3), (3.2, 12.3),
    ]},
    "Urubici": {"cor": "#B4D900", "label": (3.8, 13.2), "coords": [
        (3.3, 12.8), (3.4, 13.5), (3.7, 14), (4.0, 13.9), (4.2, 13.3), (4.3, 12.5),
        (4.0, 11.8), (3.6, 11.7), (3.4, 12.2),
    ]},
}
CAMPOS_D_ZRY_SR = {
    "Paranapanema": {"cor": "#B5B5B5", "label": (280, 5.0), "coords": [
        (150, 4.7), (160, 5.2), (200, 5.5), (280, 5.4), (350, 5.0), (400, 4.7),
        (420, 4.5), (350, 4.4), (250, 4.4), (180, 4.5),
    ]},
    "Pitanga": {"cor": "#E6951A", "label": (520, 6.8), "coords": [
        (380, 5.5), (400, 6.2), (430, 7.0), (480, 7.6), (550, 7.8), (620, 7.5),
        (650, 7.0), (630, 6.2), (560, 5.7), (480, 5.5), (420, 5.4),
    ]},
    "Urubici": {"cor": "#B4D900", "label": (820, 8.3), "coords": [
        (550, 7.6), (570, 8.2), (620, 8.7), (700, 9.0), (800, 9.1), (900, 8.8),
        (1000, 8.3), (1080, 7.8), (1100, 7.3), (1000, 7.3), (880, 7.4), (750, 7.5),
        (650, 7.4), (590, 7.4),
    ]},
}
CAMPOS_E_ZRY_TIZR = {
    "Gramado": {"cor": "#D62728", "label": (50, 4.6), "coords": [
        (33, 3.9), (35, 4.3), (40, 4.9), (46, 5.3), (52, 5.5), (58, 5.3), (62, 4.9),
        (65, 4.4), (68, 3.9), (63, 3.6), (55, 3.6), (45, 3.7), (38, 3.7),
    ]},
    "Esmeralda": {"cor": "#10962F", "label": (70, 3.7), "coords": [
        (58, 3.6), (60, 4.0), (65, 4.2), (70, 4.3), (76, 4.1), (80, 3.8), (83, 3.5),
        (80, 3.2), (74, 3.1), (66, 3.15), (60, 3.3),
    ]},
    "Ribeira": {"cor": "#8791D0", "label": (85, 3.9), "coords": [
        (63, 3.7), (65, 4.0), (70, 4.3), (78, 4.4), (88, 4.3), (96, 4.1), (100, 3.8),
        (97, 3.5), (88, 3.4), (78, 3.4), (70, 3.45),
    ]},
}
CAMPOS_F_SR_TIY = {
    "Gramado": {"cor": "#D62728", "label": (225, 240), "coords": [
        (150, 155), (155, 190), (165, 240), (180, 280), (200, 315), (225, 330),
        (250, 320), (270, 290), (285, 250), (295, 210), (300, 175), (285, 150),
        (255, 140), (220, 140), (185, 145), (165, 148),
    ]},
    "Esmeralda": {"cor": "#10962F", "label": (270, 175), "coords": [
        (220, 150), (225, 175), (240, 200), (260, 213), (285, 210), (305, 195),
        (315, 175), (310, 155), (290, 143), (260, 140), (235, 143),
    ]},
    "Ribeira": {"cor": "#8791D0", "label": (430, 260), "coords": [
        (270, 175), (280, 220), (300, 270), (330, 310), (370, 335), (420, 342),
        (470, 330), (510, 300), (535, 260), (545, 220), (535, 190), (505, 168),
        (460, 158), (410, 158), (360, 162), (315, 167),
    ]},
}


def montar_diagrama_campos(campos, titulo, x_titulo, y_titulo, x_range, y_range, suavizar=True, nota=NOTA_APROXIMADO):
    fig = go.Figure()
    campos_traces, anotacoes = montar_shapes_campos(campos, suavizar=suavizar)
    for t in campos_traces:
        fig.add_trace(t)
    fig.update_layout(annotations=anotacoes)
    fig.update_xaxes(title_text=x_titulo, range=x_range)
    fig.update_yaxes(title_text=y_titulo, range=y_range)
    return tema_grafico(fig, titulo, nota=nota)


def main():
    print("Carregando pontos de campo (Taió)...")
    registros_campo = carregar_campo()
    print(f"  {len(registros_campo)} pontos")

    # ------------------------------------------------------------
    # MAPA -- Leaflet (tiles reais, zoom/pan continuo), nao mais Plotly com
    # imagem estatica esticada -- mesma tecnica do webmap dedicado
    # (gerar_webmap_taio.py).
    # ------------------------------------------------------------
    transformer_wgs84 = Transformer.from_crs("EPSG:31982", "EPSG:4326", always_xy=True)
    campo_coord = [r for r in registros_campo if r["x"] is not None]
    if campo_coord:
        lons, lats = transformer_wgs84.transform([r["x"] for r in campo_coord], [r["y"] for r in campo_coord])
        for r, lon, lat in zip(campo_coord, lons, lats):
            r["lat"], r["lon"] = float(lat), float(lon)
    print(f"  {len(campo_coord)}/{len(registros_campo)} com coordenadas válidas pro mapa")

    lats = [r["lat"] for r in campo_coord]
    lons = [r["lon"] for r in campo_coord]
    s, n, w, e = min(lats), max(lats), min(lons), max(lons)
    folga_y, folga_x = (n - s) * 0.1, (e - w) * 0.1
    bounds_taio = [[s - folga_y, w - folga_x], [n + folga_y, e + folga_x]]

    geojson_campo = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature", "geometry": {"type": "Point", "coordinates": [r["lon"], r["lat"]]},
            "properties": {"id": r["id"], "cor": r["cor_mapa"], "popup": f"{r['nome']}<br>{r['litologia']}"},
        } for r in campo_coord],
    }

    # mapa geologico real (CPRM, 9 formacoes) -- mesma camada usada no
    # webmap dedicado (gerar_webmap_taio.py), pedido explicito do usuario
    # ("faltou o mapa geológico junto").
    formacoes_path = BASE.parent.parent / "2_Banco_de_Dados" / "saida_processada" / "formacoes_cprm_poligonos.geojson"
    geojson_formacoes = None
    if formacoes_path.exists():
        gdf_formacoes = gpd.read_file(formacoes_path).to_crs(4326)
        gdf_formacoes["popup"] = gdf_formacoes["formacao"]
        geojson_formacoes = json.loads(gdf_formacoes[["formacao", "cor", "popup", "geometry"]].to_json())
        print(f"  mapa geológico: {len(gdf_formacoes)} formações")

    # ------------------------------------------------------------
    # diagramas geoquimicos -- todos so com campos de referencia (sem
    # pontos, ver nota no topo do arquivo)
    # ------------------------------------------------------------
    fig_tas = montar_tas()
    fig_afm = montar_afm()
    fig_shand = montar_shand()
    fig_mgo_tio2 = montar_diagrama_campos(CAMPOS_MGO_TIO2, "MgO × TiO₂ (tipos de magma)", "MgO (% em peso)", "TiO₂ (% em peso)", [2, 10], [0.5, 4.8])
    fig_tiy_tizr = montar_diagrama_campos(CAMPOS_TIY_TIZR, "Ti/Y × Ti/Zr (tipos de magma)", "Ti/Y", "Ti/Zr", [100, 700], [30, 100])
    fig_fe2o3 = montar_diagrama_campos(CAMPOS_C_FE2O3_TIO2, "Fe₂O₃(t) × TiO₂ — grupo Alto-Ti", "TiO₂", "Fe₂O₃ (t)", [1, 5], [10, 19])
    fig_zry_sr = montar_diagrama_campos(CAMPOS_D_ZRY_SR, "Zr/Y × Sr — grupo Alto-Ti", "Sr", "Zr/Y", [0, 1200], [4, 10])
    fig_zry_tizr = montar_diagrama_campos(CAMPOS_E_ZRY_TIZR, "Zr/Y × Ti/Zr — grupo Baixo-Ti", "Ti/Zr", "Zr/Y", [30, 100], [3, 7])
    fig_sr_tiy = montar_diagrama_campos(CAMPOS_F_SR_TIY, "Sr × Ti/Y — grupo Baixo-Ti", "Ti/Y", "Sr", [100, 600], [100, 400])

    # ------------------------------------------------------------
    # estatisticas do catalogo de campo (Taio)
    # ------------------------------------------------------------
    contagem_lito = Counter(r["litologia"] for r in registros_campo)
    itens_lito = sorted(contagem_lito.items(), key=lambda kv: -kv[1])
    fig_lito_taio = tema_grafico(go.Figure(go.Bar(
        x=[k for k, _ in itens_lito], y=[v for _, v in itens_lito],
        marker=dict(color=[CORES_LITOLOGIA_CAMPO.get(k, COR_LITOLOGIA_PADRAO) for k, _ in itens_lito]),
    )), "Distribuição por litologia")
    fig_lito_taio.update_xaxes(tickangle=-35)

    contagem_ti_taio = Counter((r["classificacao_ti"] or "Sem dado") for r in registros_campo)
    ordem_ti = ["Alto", "Baixo", "Sem dado"]
    itens_ti = [(k, contagem_ti_taio.get(k, 0)) for k in ordem_ti if contagem_ti_taio.get(k, 0)]
    fig_ti_taio = tema_grafico(go.Figure(go.Bar(
        x=[k for k, _ in itens_ti], y=[v for _, v in itens_ti],
        marker=dict(color=["#E67E22", "#2E86C1", COR_TI_INDEFINIDO][:len(itens_ti)]),
    )), "Classificação de Ti")

    print("Montando HTML final...")
    # a biblioteca plotly.js so precisa ser embutida uma vez -- na PRIMEIRA
    # figura que aparece na pagina.
    html_tas = pio.to_html(fig_tas, full_html=False, include_plotlyjs=True, div_id="grafico-tas")
    html_afm = pio.to_html(fig_afm, full_html=False, include_plotlyjs=False, div_id="grafico-afm")
    html_shand = pio.to_html(fig_shand, full_html=False, include_plotlyjs=False, div_id="grafico-shand")
    html_mgo_tio2 = pio.to_html(fig_mgo_tio2, full_html=False, include_plotlyjs=False, div_id="grafico-mgo-tio2")
    html_tiy_tizr = pio.to_html(fig_tiy_tizr, full_html=False, include_plotlyjs=False, div_id="grafico-tiy-tizr")
    html_fe2o3 = pio.to_html(fig_fe2o3, full_html=False, include_plotlyjs=False, div_id="grafico-fe2o3")
    html_zry_sr = pio.to_html(fig_zry_sr, full_html=False, include_plotlyjs=False, div_id="grafico-zry-sr")
    html_zry_tizr = pio.to_html(fig_zry_tizr, full_html=False, include_plotlyjs=False, div_id="grafico-zry-tizr")
    html_sr_tiy = pio.to_html(fig_sr_tiy, full_html=False, include_plotlyjs=False, div_id="grafico-sr-tiy")
    html_lito_taio = pio.to_html(fig_lito_taio, full_html=False, include_plotlyjs=False, div_id="grafico-lito-taio")
    html_ti_taio = pio.to_html(fig_ti_taio, full_html=False, include_plotlyjs=False, div_id="grafico-ti-taio")

    dados_js = json.dumps(registros_campo, ensure_ascii=False)
    logo_b64 = logo_base64()

    linhas_tabela = []
    for r in registros_campo:
        linhas_tabela.append(f"""
        <tr class="linha-dado" data-id="{r['id']}"
            data-busca="{(r['nome'] + ' ' + (r['litologia'] or '')).lower()}">
            <td><span class="nome-ponto">{r['nome']}</span></td>
            <td>{r['litologia'] or '—'}</td>
            <td>{r['classificacao_ti'] or '—'}</td>
        </tr>""")
    tabela_html = "".join(linhas_tabela)

    html_final = f"""<!DOCTYPE html>
<html lang="pt-br">
<head>
<meta charset="utf-8">
<link rel="icon" type="image/png" href="assets/favicon.png">
<link rel="shortcut icon" href="assets/favicon.ico">
<link rel="apple-touch-icon" href="assets/apple-touch-icon.png">
<title>Dashboard Geoquímico — Taió</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"
      integrity="sha256-p4NxAoJBhIIN+hmNHrzRCf9tD/miZyoHS5obTRR9BMY=" crossorigin="">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"
        integrity="sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo=" crossorigin=""></script>
<style>
  * {{ box-sizing: border-box; }}
  :root {{
    --bg: {MARCA_NAVY}; --painel: {COR_PAINEL}; --texto: {MARCA_CINZA_CLARO};
    --borda-fraca: #3a3f52; --input-bg: {MARCA_NAVY};
  }}
  body.tema-claro {{
    --bg: {MARCA_CINZA_CLARO}; --painel: #FFFFFF; --texto: {MARCA_NAVY};
    --borda-fraca: #D8D8E2; --input-bg: #FFFFFF;
  }}
  body {{ margin: 0; background: var(--bg); color: var(--texto); font-family: {MARCA_FONTE}; transition: background 0.2s, color 0.2s; }}
  header {{ display: flex; align-items: center; justify-content: space-between; padding: 14px 24px; border-bottom: 1px solid {MARCA_ROXO}; gap: 16px; }}
  header h1 {{ font-size: 20px; margin: 0; }}
  header h1 b {{ color: {MARCA_ROXO}; }}
  header .lado-direito {{ display: flex; align-items: center; gap: 12px; }}
  header img {{ width: 52px; height: 52px; border-radius: 50%; border: 2px solid {MARCA_ROXO}; box-shadow: 0 0 10px rgba(123,47,255,0.6); }}
  .btn-tema {{
    background: {MARCA_ROXO_ESCURO}; color: {MARCA_CINZA_CLARO}; border: 1.5px solid {MARCA_ROXO};
    border-radius: 6px; padding: 7px 12px; font-family: {MARCA_FONTE}; font-size: 12px; cursor: pointer;
  }}
  body.tema-claro .btn-tema {{ background: #EDE3FF; color: {MARCA_NAVY}; }}
  .btn-tema.ativo {{ opacity: 1; font-weight: 700; }}
  .btn-tema:not(.ativo) {{ opacity: 0.55; }}
  .layout {{ display: grid; grid-template-columns: 300px 1fr 1fr; gap: 12px; padding: 12px; height: calc(100vh - 90px); min-height: 600px; }}
  .painel {{ background: var(--painel); border: 1px solid {MARCA_ROXO}; border-radius: 8px; overflow: hidden; display: flex; flex-direction: column; }}
  .painel h2 {{ font-size: 13px; margin: 0; padding: 10px 14px; border-bottom: 1px solid var(--borda-fraca); color: var(--texto); opacity: 0.85; text-transform: uppercase; letter-spacing: 0.05em; }}
  #busca {{ margin: 10px 14px 0 14px; padding: 7px 10px; border-radius: 6px; border: 1px solid {MARCA_ROXO}; background: var(--input-bg); color: var(--texto); font-family: {MARCA_FONTE}; }}
  .contagem {{ font-size: 11px; opacity: 0.6; padding: 6px 14px 0 14px; }}
  .lista-scroll {{ overflow-y: auto; flex: 1; padding: 0 8px 8px 8px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
  th {{ position: sticky; top: 0; background: var(--painel); text-align: left; padding: 6px 6px; border-bottom: 1px solid {MARCA_ROXO}; opacity: 0.75; font-weight: 500; }}
  td {{ padding: 6px 6px; border-bottom: 1px solid var(--borda-fraca); }}
  .linha-dado {{ cursor: pointer; }}
  .linha-dado:hover {{ background: rgba(123,47,255,0.15); }}
  .linha-dado.selecionada {{ background: rgba(123,47,255,0.35); }}
  .nome-ponto {{
    display: inline-block; padding: 2px 7px; border-radius: 5px;
    border: 1px solid rgba(123,47,255,0.35); background: rgba(123,47,255,0.08);
  }}
  .linha-dado.selecionada .nome-ponto {{
    color: {MARCA_ROXO}; font-weight: 700; border-color: {MARCA_ROXO}; background: rgba(123,47,255,0.18);
  }}
  body.tema-claro .linha-dado.selecionada .nome-ponto {{ color: #5A1FBF; border-color: #5A1FBF; }}
  footer {{ text-align: center; padding: 8px; opacity: 0.55; font-size: 11px; }}
  #mapa-leaflet {{ flex: 1; }}
  .col-graficos {{ overflow-y: auto; padding: 10px; display: flex; flex-direction: column; gap: 12px; }}
  .grafico-card {{ background: var(--painel); border: 1px solid var(--borda-fraca); border-left: 4px solid {MARCA_ROXO}; border-radius: 8px; flex-shrink: 0; }}
  .leaflet-popup-content-wrapper {{ background: {MARCA_ROXO_ESCURO}; color: {MARCA_CINZA_CLARO}; border: 1px solid {MARCA_ROXO}; }}
  .leaflet-popup-tip {{ background: {MARCA_ROXO_ESCURO}; }}
  .leaflet-popup-content {{ font-family: {MARCA_FONTE}; font-size: 12px; }}
  .leaflet-control-layers {{ background: {MARCA_ROXO_ESCURO} !important; color: {MARCA_CINZA_CLARO}; border: 1px solid {MARCA_ROXO} !important; font-size: 12px; }}
  .leaflet-control-layers-toggle {{ filter: invert(1); }}
  .leaflet-bar a {{ background: {MARCA_ROXO_ESCURO}; color: {MARCA_CINZA_CLARO}; border-bottom-color: {MARCA_ROXO} !important; }}
  .leaflet-bar a:hover {{ background: {MARCA_ROXO}; }}
  #popup-info {{
    position: fixed; z-index: 2000; display: none; max-width: 320px;
    background: {MARCA_ROXO_ESCURO}; border: 1px solid {MARCA_ROXO}; border-radius: 8px;
    padding: 10px 12px; font-size: 12px; box-shadow: 0 4px 16px rgba(0,0,0,0.5); pointer-events: none;
  }}
  #popup-info b {{ color: {MARCA_CINZA_CLARO}; }}
  #popup-info .linha-popup {{ margin: 2px 0; opacity: 0.9; }}
  #popup-info .titulo-popup {{ font-size: 13px; font-weight: 700; color: {MARCA_ROXO}; margin-bottom: 4px; }}
</style>
</head>
<body>
<header>
  <h1><b>Dashboard Geoquímico</b> — Taió</h1>
  <div class="lado-direito">
    <button class="btn-tema ativo" id="btn-tema-escuro" onclick="aplicarTemaGeral('escuro')">Tema: Escuro</button>
    <button class="btn-tema" id="btn-tema-claro" onclick="aplicarTemaGeral('claro')">Tema: Claro</button>
    {f'<img src="data:image/jpeg;base64,{logo_b64}">' if logo_b64 else ''}
  </div>
</header>
<div class="layout">
  <div class="painel">
    <h2>Lista de dados</h2>
    <input id="busca" type="text" placeholder="Buscar por nome, litologia...">
    <div class="contagem" id="contagem"></div>
    <div class="lista-scroll">
      <table>
        <thead><tr><th>Nome</th><th>Litologia</th><th>Ti</th></tr></thead>
        <tbody id="corpo-tabela">{tabela_html}</tbody>
      </table>
    </div>
  </div>
  <div class="painel col-graficos">
    <div class="grafico-card">{html_tas}</div>
    <div class="grafico-card">{html_afm}</div>
    <div class="grafico-card">{html_shand}</div>
    <div class="grafico-card">{html_mgo_tio2}</div>
    <div class="grafico-card">{html_tiy_tizr}</div>
    <div class="grafico-card">{html_fe2o3}</div>
    <div class="grafico-card">{html_zry_sr}</div>
    <div class="grafico-card">{html_zry_tizr}</div>
    <div class="grafico-card">{html_sr_tiy}</div>
    <div class="grafico-card">{html_lito_taio}</div>
    <div class="grafico-card">{html_ti_taio}</div>
  </div>
  <div class="painel">
    <h2>Mapa</h2>
    <div id="mapa-leaflet"></div>
  </div>
</div>
<div id="popup-info"></div>
<footer>Criado por Afonso Henrique de Jesus</footer>
<script>
(function() {{
    var DADOS = {dados_js};
    var DADOS_POR_ID = {{}};
    DADOS.forEach(function(r) {{ DADOS_POR_ID[r.id] = r; }});

    var TODOS_GD = [
        document.getElementById('grafico-tas'), document.getElementById('grafico-afm'),
        document.getElementById('grafico-shand'), document.getElementById('grafico-mgo-tio2'),
        document.getElementById('grafico-tiy-tizr'), document.getElementById('grafico-fe2o3'),
        document.getElementById('grafico-zry-sr'), document.getElementById('grafico-zry-tizr'),
        document.getElementById('grafico-sr-tiy'),
        document.getElementById('grafico-lito-taio'), document.getElementById('grafico-ti-taio'),
    ];

    // tema claro/escuro -- moldura (fundo/eixos/legenda/botoes) muda, cores
    // dos dados (litologia, campos de literatura) ficam fixas.
    var TEMA = {{
        escuro: {{
            paper: '{MARCA_NAVY}', painel: '{COR_PAINEL}', texto: '{MARCA_CINZA_CLARO}',
            grid: '#3a3f52', zerogrid: '#3a3f52', legendBg: 'rgba(45,10,74,0.75)', botaoBg: '{MARCA_ROXO_ESCURO}',
        }},
        claro: {{
            paper: '{MARCA_CINZA_CLARO}', painel: '#FFFFFF', texto: '{MARCA_NAVY}',
            grid: '#D8D8E2', zerogrid: '#C4C4D0', legendBg: 'rgba(255,255,255,0.85)', botaoBg: '#EDE3FF',
        }},
    }};

    function temaPlotly(gd, nome) {{
        if (!gd || !gd.layout) return;
        var t = TEMA[nome];
        var patch = {{
            paper_bgcolor: t.paper, plot_bgcolor: t.painel, 'font.color': t.texto,
            'xaxis.color': t.texto, 'xaxis.gridcolor': t.grid, 'xaxis.zerolinecolor': t.zerogrid,
            'yaxis.color': t.texto, 'yaxis.gridcolor': t.grid, 'yaxis.zerolinecolor': t.zerogrid,
            'legend.bgcolor': t.legendBg, 'legend.font.color': t.texto, 'title.font.color': t.texto,
        }};
        if (gd.layout.ternary) {{
            patch['ternary.bgcolor'] = t.painel;
            patch['ternary.aaxis.color'] = t.texto; patch['ternary.aaxis.gridcolor'] = t.grid;
            patch['ternary.baxis.color'] = t.texto; patch['ternary.baxis.gridcolor'] = t.grid;
            patch['ternary.caxis.color'] = t.texto; patch['ternary.caxis.gridcolor'] = t.grid;
        }}
        Plotly.relayout(gd, patch);
    }}

    window.aplicarTemaGeral = function(nome) {{
        document.body.classList.toggle('tema-claro', nome === 'claro');
        document.getElementById('btn-tema-escuro').classList.toggle('ativo', nome === 'escuro');
        document.getElementById('btn-tema-claro').classList.toggle('ativo', nome === 'claro');
        TODOS_GD.forEach(function(gd) {{ temaPlotly(gd, nome); }});
    }};

    // ---- mapa Leaflet ----
    var mapa = L.map('mapa-leaflet', {{ zoomControl: true }});
    var satelite = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{{z}}/{{y}}/{{x}}', {{
        attribution: 'Esri World Imagery', maxZoom: 19,
    }});
    var rico = L.tileLayer('https://{{s}}.basemaps.cartocdn.com/rastertiles/voyager/{{z}}/{{x}}/{{y}}{{r}}.png', {{
        attribution: '&copy; OpenStreetMap &copy; CARTO', maxZoom: 20, subdomains: 'abcd',
    }});
    var escuro = L.tileLayer('https://{{s}}.basemaps.cartocdn.com/dark_all/{{z}}/{{x}}/{{y}}{{r}}.png', {{
        attribution: '&copy; OpenStreetMap &copy; CARTO', maxZoom: 20, subdomains: 'abcd',
    }});
    var osmPadrao = L.tileLayer('https://tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png', {{
        attribution: '&copy; OpenStreetMap contributors', maxZoom: 19,
    }});
    var relevo = L.tileLayer('https://{{s}}.tile.opentopomap.org/{{z}}/{{x}}/{{y}}.png', {{
        attribution: '&copy; OpenStreetMap contributors, SRTM &copy; OpenTopoMap (CC-BY-SA)',
        maxZoom: 17, subdomains: 'abc',
    }});
    escuro.addTo(mapa);
    mapa.fitBounds({json.dumps(bounds_taio)});

    function pontoEstilo(cor, raio) {{
        return {{ radius: raio, fillColor: cor, color: '{MARCA_CINZA_CLARO}', weight: 1, fillOpacity: 0.9 }};
    }}
    var campoLayer = L.geoJSON({json.dumps(geojson_campo, ensure_ascii=False)}, {{
        pointToLayer: function(f, latlng) {{ return L.circleMarker(latlng, pontoEstilo(f.properties.cor, 6)); }},
        onEachFeature: function(f, layer) {{
            layer.bindPopup(f.properties.popup);
            layer.on('click', function() {{ selecionarPorId(f.properties.id, 'mapa'); }});
        }},
    }}).addTo(mapa);
    var overlaysMapa = {{ "Campo (Taió)": campoLayer }};"""
    if geojson_formacoes is not None:
        html_final += f"""
    var formacoesLayer = L.geoJSON({json.dumps(geojson_formacoes, ensure_ascii=False)}, {{
        style: function(f) {{ return {{ color: '#000', weight: 0.5, fillColor: f.properties.cor, fillOpacity: 0.5 }}; }},
        onEachFeature: function(f, layer) {{ layer.bindPopup(f.properties.popup); }},
    }}).addTo(mapa);
    overlaysMapa["Mapa geológico real (CPRM)"] = formacoesLayer;"""
    html_final += f"""
    L.control.layers(
        {{ "Escuro (CartoDB Dark)": escuro, "Rico (CartoDB Voyager)": rico, "Satélite (Esri)": satelite, "Relevo/Topográfico": relevo, "OSM Padrão": osmPadrao }},
        overlaysMapa, {{ collapsed: false }}
    ).addTo(mapa);
    L.control.scale({{ metric: true, imperial: false }}).addTo(mapa);

    var destaqueMapa = L.circleMarker([0, 0], {{
        radius: 14, color: '{MARCA_ROXO}', weight: 3, fillOpacity: 0, opacity: 0,
    }}).addTo(mapa);

    function selecionarPorId(id, origemClique) {{
        document.querySelectorAll('.linha-dado.selecionada').forEach(function(el) {{ el.classList.remove('selecionada'); }});
        var linha = document.querySelector('.linha-dado[data-id="' + id + '"]');
        if (linha) {{
            linha.classList.add('selecionada');
            if (origemClique !== 'lista') linha.scrollIntoView({{block: 'nearest'}});
        }}
        var r = DADOS_POR_ID[id];
        if (!r) return;
        if (r.lat !== undefined && r.lon !== undefined) {{
            destaqueMapa.setLatLng([r.lat, r.lon]);
            destaqueMapa.setStyle({{ opacity: 1 }});
            if (origemClique !== 'mapa') {{
                // da um leve zoom junto (nao so pan) -- se ja estiver bem
                // perto (zoom alto), nao forca zoom pra fora, so aproxima
                // quando fizer sentido (zoom atual menor que o alvo).
                var zoomAlvo = Math.max(mapa.getZoom(), 14);
                mapa.flyTo([r.lat, r.lon], zoomAlvo, {{ duration: 0.6 }});
            }}
        }} else {{
            destaqueMapa.setStyle({{ opacity: 0 }});
        }}
    }}

    document.getElementById('corpo-tabela').addEventListener('click', function(ev) {{
        var linha = ev.target.closest('.linha-dado');
        if (!linha) return;
        selecionarPorId(linha.getAttribute('data-id'), 'lista');
    }});

    // popup com os detalhes do ponto ao passar o mouse na linha da lista
    var popup = document.getElementById('popup-info');
    function montarPopup(r) {{
        var linhas = [];
        linhas.push('<div class="titulo-popup">' + r.nome + '</div>');
        if (r.litologia) linhas.push('<div class="linha-popup"><b>Litologia:</b> ' + r.litologia + '</div>');
        if (r.classificacao_ti) linhas.push('<div class="linha-popup"><b>Ti:</b> ' + r.classificacao_ti + '</div>');
        if (r.tipo_ponto) linhas.push('<div class="linha-popup"><b>Tipo:</b> ' + r.tipo_ponto + '</div>');
        if (r.qualidade) linhas.push('<div class="linha-popup"><b>Qualidade:</b> ' + r.qualidade + '</div>');
        if (r.x !== null && r.x !== undefined) linhas.push('<div class="linha-popup"><b>UTM:</b> ' + Math.round(r.x) + ', ' + Math.round(r.y) + '</div>');
        if (r.descricao) linhas.push('<div class="linha-popup" style="margin-top:4px; opacity:0.75;">' + r.descricao + '</div>');
        return linhas.join('');
    }}
    document.getElementById('corpo-tabela').addEventListener('mouseover', function(ev) {{
        var linha = ev.target.closest('.linha-dado');
        if (!linha) return;
        var r = DADOS_POR_ID[linha.getAttribute('data-id')];
        if (!r) return;
        popup.innerHTML = montarPopup(r);
        popup.style.display = 'block';
    }});
    document.getElementById('corpo-tabela').addEventListener('mousemove', function(ev) {{
        if (popup.style.display !== 'block') return;
        var x = ev.clientX + 16, y = ev.clientY + 12;
        if (x + 330 > window.innerWidth) x = ev.clientX - 336;
        popup.style.left = x + 'px';
        popup.style.top = y + 'px';
    }});
    document.getElementById('corpo-tabela').addEventListener('mouseout', function(ev) {{
        var linha = ev.target.closest('.linha-dado');
        if (!linha) return;
        popup.style.display = 'none';
    }});

    // filtro de busca
    function aplicarFiltro() {{
        var termo = document.getElementById('busca').value.toLowerCase();
        var visiveis = 0;
        document.querySelectorAll('.linha-dado').forEach(function(el) {{
            var mostra = !termo || el.getAttribute('data-busca').indexOf(termo) !== -1;
            el.style.display = mostra ? '' : 'none';
            if (mostra) visiveis++;
        }});
        document.getElementById('contagem').textContent = visiveis + ' de ' + DADOS.length + ' registros';
    }}
    document.getElementById('busca').addEventListener('input', aplicarFiltro);
    aplicarFiltro();
}})();
</script>
</body>
</html>
"""
    OUT_HTML.write_text(html_final, encoding="utf-8")
    print(f"Salvo em: {OUT_HTML}")


if __name__ == "__main__":
    main()
