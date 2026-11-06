"""Dashboard geoquimico -- terceiro produto web (lista + mapa + diagramas),
100% escopo Taio (a versao anterior misturava com uma planilha de geoquimica
da PMP/Florianopolis -- removida por pedido explicito do usuario; o produto
PMP vai virar um dashboard SEPARADO, sem afetar este aqui).

Dados: 308 pontos de campo de `pontos_unificados_completo.gpkg` (mesmo
catalogo unificado usado nos outros 3 produtos) + 41 amostras de geoquimica
bruta resgatadas em `2_Banco_de_Dados/QMC_TAIO_TODOS` (26 de
afloramento/campo, cruzadas com o ponto de campo por SAMPLE NAME =
id_original, + 15 de furo na soleira Bela Vista + 3 amostras de referencia
com nomenclatura ainda a confirmar) -- ver
`QMC_TAIO_TODOS/Geoquimica_Taio_Cruzada.xlsx` pro cruzamento/analise previa
completa.

Os diagramas geoquimicos (TAS, AFM, Shand's Index, e 6 campos de magma-tipo
digitalizados de figura da literatura) agora mostram as 41 amostras reais
plotadas por cima dos campos de referencia (que continuam aproximados,
digitalizados a olho -- ver nota em cada figura). Cores: laranja = Alto-Ti,
azul = Baixo-Ti (mesmo criterio TiO2 >= 2,0% em peso ja usado no
pontos_unificados.csv, validado sem divergencia contra as 21 amostras em
comum). O mapa ganhou uma camada extra "Geoquímica (amostras brutas)" com a
localização de cada amostra (a do furo repete a mesma boca de sondagem, só
a profundidade muda).

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
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from shapely.geometry import Point
from shapely.prepared import prep

BASE = Path(__file__).parent
PONTOS_CAMPO_GPKG = (
    BASE.parent.parent / "2_Banco_de_Dados" / "Unificação" / "GPKG_Novos" / "pontos_unificados_completo.gpkg"
)
GEOQUIMICA_CSV = BASE.parent.parent / "2_Banco_de_Dados" / "QMC_TAIO_TODOS" / "geoquimica_dashboard.csv"
# litologia_processada.shp (ETL: 2_Banco_de_Dados/scripts_etl/processar_litologia_atualizada.py)
# substitui o mapa geologico real (CPRM) + o poligon_intrusiva.shp antigo --
# um shp so, com sill/dique redigitalizados (coluna "formacao") e as 6
# formacoes sedimentares, incl. deposito quaternario (coluna "tipo" ==
# "sedimentar"/"intrusiva").
LITOLOGIA_ATUALIZADA = BASE.parent.parent / "2_Banco_de_Dados" / "dados_base" / "litologia_processada.shp"
TOPO_NPY = BASE.parent / "dados_entrada" / "topografia_drone" / "topografia_xyz.npy"

# mesmo plano de mergulho regional + erosao contra o relevo real ja usado pra
# desenhar a camada Serra Alta e a base do sill no visualizador 3D/secao 2D
# (ver elevacao_serra_alta em gerar_visualizador_3d.py) -- reaproveitado aqui
# so pra estimar o VOLUME do sill (topo real menos topo da Serra Alta).
TREND_A, TREND_B = 0.01034, -0.00025
TREND_X0, TREND_Y0 = 592300.0, 7015058.8
Z_REF_TILT = 1053.5
PROFUNDIDADE_SERRA_ALTA = 350.0
ESPESSURA_MINIMA_SILL = 27.0  # mediana real de campo -- mesmo piso do visualizador 3D
ESPESSURA_MEDIANA_DIQUE = 6.0  # mediana real de campo (05_gerar_solidos_visualizacao.py) -- o
# dique nao tem uma base geologica equivalente a Serra Alta do sill (corta a sequencia toda
# quase na vertical), entao o volume usa essa espessura tabular constante x area
LOGO_PATH = BASE / "assets" / "logo_gstech.jpg"
OUT_HTML = BASE / "dashboard_geoquimico.html"

COR_TI_ALTO = "#E67E22"
COR_TI_BAIXO = "#2E86C1"

# identidade visual GS Tech -- mesma paleta dos outros produtos, manter em sincronia.
MARCA_ROXO_ESCURO = "#2D0A4A"
MARCA_ROXO = "#7B2FFF"
MARCA_NAVY = "#1B1F2E"
MARCA_CINZA_CLARO = "#F2F2F2"
MARCA_FONTE = "Montserrat, Arial, sans-serif"
COR_PAINEL = "#262B3D"

COR_SILL = "#49B18C"
COR_DIQUE = "#AFEC7C"
NOMES_CAMADAS = ["Teresina", "Serra Alta", "Irati", "Palermo", "Rio Bonito"]
CORES_CAMADAS = ["#F8B6A0", "#FB8D74", "#C8625D", "#FF7757", "#E3644F"]  # paleta do mapa geologico (Teresina, Serra Alta, Irati, Palermo, Rio Bonito)
COR_QUATERNARIO = "#FFFBC2"
CORES_LITOLOGIA_CAMPO = {
    "sill_diabasio": COR_SILL, "sill_diabasio_cprm": COR_SILL,
    "dique": COR_DIQUE, "dique_cprm": COR_DIQUE,
    "encaixante_teresina": "#F8B6A0", "encaixante_serra_alta": "#FB8D74",
    "encaixante_irati": "#C8625D", "encaixante_palermo": "#FF7757",
    "encaixante_rio_bonito": "#E3644F", "encaixante_sedimentar": "#999999",
}
COR_LITOLOGIA_PADRAO = "#999999"
COR_TI_INDEFINIDO = "#9A9A9A"


NOME_PONTO = {}  # ponto_id -> nome padrao ITC (preenchido em main)


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
            "ponto_id": row.ponto_id,
            "nome": row.nome_itc if "nome_itc" in gdf.columns and pd.notna(row.nome_itc) else row.ponto_id,  # padrao ITC
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
# 1b. geoquimica bruta resgatada (QMC_TAIO_TODOS) -- 41 amostras (26 de
#     afloramento/campo + 15 de furo na soleira + 3 amostras de referencia
#     de nomenclatura ainda a confirmar), cruzadas com pontos_unificados
#     por SAMPLE NAME = id_original. Ver
#     2_Banco_de_Dados/QMC_TAIO_TODOS/Geoquimica_Taio_Cruzada.xlsx pro
#     detalhamento/analise previa completa.
# ======================================================================
PESO_MOLAR = {"Al2O3": 101.96, "CaO": 56.08, "Na2O": 61.98, "K2O": 94.20}


def carregar_geoquimica():
    if not GEOQUIMICA_CSV.exists():
        return []
    df = pd.read_csv(GEOQUIMICA_CSV)
    registros = []
    for row in df.itertuples():
        cor = COR_TI_ALTO if row.classificacao_ti == "Alto-Ti" else COR_TI_BAIXO
        # moles (oxido/peso molar) pro indice de Shand -- A/CNK e A/NK
        mol_al = row.Al2O3 / PESO_MOLAR["Al2O3"]
        mol_ca = row.CaO / PESO_MOLAR["CaO"]
        mol_na = row.Na2O / PESO_MOLAR["Na2O"]
        mol_k = row.K2O / PESO_MOLAR["K2O"]
        a_cnk = mol_al / (mol_ca + mol_na + mol_k) if (mol_ca + mol_na + mol_k) else None
        a_nk = mol_al / (mol_na + mol_k) if (mol_na + mol_k) else None
        # AFM (wt%, normalizado pra somar 100 -- convencao usual do diagrama)
        a_wt, f_wt, m_wt = row.Na2O + row.K2O, row.Fe2O3T, row.MgO
        soma_afm = a_wt + f_wt + m_wt
        ponto_id = row.ponto_id if pd.notna(row.ponto_id) else None
        registros.append({
            "id": f"GEOQ-{row.amostra}",
            # id_display: pra amostra ja associada a um ponto de campo (tem
            # ponto_id), clique no mapa/grafico seleciona a linha do PONTO DE
            # CAMPO na lista (mesma amostra, sem duplicar linha) -- so as
            # amostras sem ponto de campo (furo/referencia) usam o proprio id.
            "id_display": f"CAMPO-{ponto_id}" if ponto_id else f"GEOQ-{row.amostra}",
            "amostra": row.amostra,
            "origem": row.origem,
            "classificacao_ti": row.classificacao_ti,
            "cor": cor,
            "ponto_id": ponto_id,
            "litologia": row.litologia if pd.notna(row.litologia) else None,
            "lat": float(row.lat) if pd.notna(row.lat) else None,
            "lon": float(row.lon) if pd.notna(row.lon) else None,
            "utm_e": float(row.utm_e) if pd.notna(row.utm_e) else None,
            "utm_n": float(row.utm_n) if pd.notna(row.utm_n) else None,
            "profundidade_m": float(row.profundidade_m) if pd.notna(row.profundidade_m) else None,
            "SiO2": row.SiO2, "TiO2": row.TiO2, "Al2O3": row.Al2O3, "FeOT": row.FeOT, "Fe2O3T": row.Fe2O3T,
            "CaO": row.CaO, "MgO": row.MgO, "K2O": row.K2O, "Na2O": row.Na2O, "P2O5": row.P2O5,
            "Sr": row.Sr, "Zr": row.Zr, "Y": row.Y, "Ti_Zr": row.Ti_Zr, "Zr_Y": row.Zr_Y,
            "Ti_Y": row.Ti_Zr * row.Zr_Y if pd.notna(row.Ti_Zr) and pd.notna(row.Zr_Y) else None,
            "tas_x": row.SiO2, "tas_y": row.Na2O + row.K2O,
            # convencao do diagrama AFM ja usado neste arquivo: a=F(topo), b=A(esquerda), c=M(direita)
            "afm_a_F": (f_wt / soma_afm * 100) if soma_afm else None,
            "afm_b_A": (a_wt / soma_afm * 100) if soma_afm else None,
            "afm_c_M": (m_wt / soma_afm * 100) if soma_afm else None,
            "shand_acnk": a_cnk, "shand_ank": a_nk,
        })
    return registros


# ======================================================================
# 1c. resumo dos corpos intrusivos (area digitalizada em planta +
#     classificacao Ti das amostras) -- pro quadro estatico embaixo da lista.
# ======================================================================
def elevacao_serra_alta(x, y, elevacao_fn):
    """Cota do topo da Serra Alta num ponto (x,y) -- mesma formula usada em
    gerar_visualizador_3d.py, reaproveitada so pra estimar volume aqui."""
    plano = Z_REF_TILT + TREND_A * (x - TREND_X0) + TREND_B * (y - TREND_Y0) - PROFUNDIDADE_SERRA_ALTA
    return min(plano, elevacao_fn(x, y))


def calcular_volume_sill_m3(grupo_sill, elevacao_fn, passo=80.0):
    """Volume do sill (m3) integrando numericamente (topo real - base
    geologica) numa grade regular sobre cada lobo -- mesma base "topo da
    Serra Alta com piso minimo de 27m" usada no visualizador 3D/secao 2D
    (ver elevacao_serra_alta), so que aqui vira um numero agregado em vez de
    uma malha 3D."""
    volume_total = 0.0
    celula_m2 = passo * passo
    for poligono in grupo_sill.geometry:
        partes = poligono.geoms if poligono.geom_type == "MultiPolygon" else [poligono]
        for parte in partes:
            minx, miny, maxx, maxy = parte.bounds
            xs = np.arange(minx + passo / 2, maxx, passo)
            ys = np.arange(miny + passo / 2, maxy, passo)
            if len(xs) == 0 or len(ys) == 0:
                xs, ys = np.array([(minx + maxx) / 2]), np.array([(miny + maxy) / 2])
            mascara_pt = prep(parte)
            for x in xs:
                for y in ys:
                    if not mascara_pt.contains(Point(x, y)):
                        continue
                    z_topo = elevacao_fn(x, y)
                    z_base = min(elevacao_serra_alta(x, y, elevacao_fn), z_topo - ESPESSURA_MINIMA_SILL)
                    volume_total += (z_topo - z_base) * celula_m2
    return volume_total


PROFUNDIDADES_CAMADAS = [0.0, 350.0, 430.0, 485.0, 585.0, 854.0]  # topo/base de Teresina..Rio Bonito
QUATERNARIO_LIMIAR_REAL = 450.0  # mesmo limiar/espessura do visualizador 3D e da secao 2D
QUATERNARIO_ESPESSURA_REAL = 30.0


def _cor_texto_sobre(hex_cor):
    r, g, b = (int(hex_cor.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4))
    return "#1B1F2E" if (0.299 * r + 0.587 * g + 0.114 * b) > 140 else "#F2F2F2"


def calcular_colunas(registros):
    """Coluna estratigrafica ESTIMADA em cada ponto (como um furo virtual):
    mesmo modelo da secao 2D/visualizador 3D -- contato[0] = relevo real,
    contatos de subsuperficie = plano de mergulho regional menos a
    profundidade acumulada, "capados" pelo relevo (erosao); quaternario so
    onde o relevo fica abaixo do limiar; sill (poligono) substitui o trecho
    superior da coluna ate a base calculada; dique (poligono) vira uma barra
    lateral atravessando a coluna toda. Devolve {id: {...}} pro JS."""
    if not TOPO_NPY.exists() or not LITOLOGIA_ATUALIZADA.exists():
        return {}, 0.0
    xyz = np.load(TOPO_NPY)
    linear = LinearNDInterpolator(xyz[:, :2], xyz[:, 2])
    nearest = NearestNDInterpolator(xyz[:, :2], xyz[:, 2])

    def elevacao_fn(x, y):
        z = linear(x, y)
        return float(z) if not np.isnan(z) else float(nearest(x, y))

    gdf = gpd.read_file(LITOLOGIA_ATUALIZADA)
    sill_geom = gdf[gdf["formacao"] == "Soleira"].geometry.union_all()
    dique_geom = gdf[gdf["formacao"] == "Dique"].geometry.union_all()
    cores = dict(zip(NOMES_CAMADAS, CORES_CAMADAS))

    colunas, max_total = {}, 0.0
    for r in registros:
        if r.get("x") is None or r.get("y") is None:
            continue
        x, y = float(r["x"]), float(r["y"])
        z = elevacao_fn(x, y)
        tilt = Z_REF_TILT + TREND_A * (x - TREND_X0) + TREND_B * (y - TREND_Y0)
        contatos = [z if m == 0 else min(tilt - PROFUNDIDADES_CAMADAS[m], z) for m in range(len(PROFUNDIDADES_CAMADAS))]
        ponto = Point(x, y)
        no_sill = sill_geom.contains(ponto)
        no_dique = dique_geom.contains(ponto)

        segs = []  # [nome, topo_m, base_m, cor]
        if no_sill:
            base_sill = max(ESPESSURA_MINIMA_SILL, z - elevacao_serra_alta(x, y, elevacao_fn))
            segs.append(["Soleira (diabásio)", 0.0, base_sill, COR_SILL])
            piso = base_sill
        elif z <= QUATERNARIO_LIMIAR_REAL:
            segs.append(["Depósito quaternário", 0.0, QUATERNARIO_ESPESSURA_REAL, COR_QUATERNARIO])
            piso = QUATERNARIO_ESPESSURA_REAL
        else:
            piso = 0.0
        for k, nome in enumerate(NOMES_CAMADAS):
            topo, base = z - contatos[k], z - contatos[k + 1]
            topo = max(topo, piso)
            if base - topo > 0.5:
                segs.append([nome, topo, base, cores[nome]])
        total = z - contatos[-1]
        max_total = max(max_total, total)
        colunas[r["id"]] = {
            "z": round(z), "total": round(total),
            "segs": [[n, round(t), round(b), c, _cor_texto_sobre(c)] for n, t, b, c in segs],
            "dique": bool(no_dique),
            "d_sill": None if no_sill else round(sill_geom.distance(ponto)),
            "d_dique": None if no_dique else round(dique_geom.distance(ponto)),
        }
    return colunas, max_total


def carregar_resumo_corpos(registros_geoq):
    resumo = {"sill_diabasio": {"nome": "Soleira (sill)", "area_km2": 0.0, "n_poligonos": 0, "alto": 0, "baixo": 0, "volume_km3": None},
              "dique": {"nome": "Dique", "area_km2": 0.0, "n_poligonos": 0, "alto": 0, "baixo": 0, "volume_km3": None}}
    if LITOLOGIA_ATUALIZADA.exists():
        gdf = gpd.read_file(LITOLOGIA_ATUALIZADA)
        mapa_tipo = {"Soleira": "sill_diabasio", "Dique": "dique"}
        for tipo, grupo in gdf.groupby("formacao"):
            chave = mapa_tipo.get(tipo)
            if not chave:
                continue
            area_m2 = float(grupo.geometry.area.sum())
            resumo[chave]["area_km2"] = area_m2 / 1e6
            resumo[chave]["n_poligonos"] = int(len(grupo))
            if chave == "dique":
                resumo[chave]["volume_km3"] = area_m2 * ESPESSURA_MEDIANA_DIQUE / 1e9
            elif chave == "sill_diabasio" and TOPO_NPY.exists():
                xyz = np.load(TOPO_NPY)
                linear = LinearNDInterpolator(xyz[:, :2], xyz[:, 2])
                nearest = NearestNDInterpolator(xyz[:, :2], xyz[:, 2])

                def elevacao_fn(x, y):
                    z = linear(x, y)
                    return float(z) if not np.isnan(z) else float(nearest(x, y))

                resumo[chave]["volume_km3"] = calcular_volume_sill_m3(grupo, elevacao_fn) / 1e9
    # amostras de furo (Soleira Bela Vista) nao tem litologia preenchida no
    # CSV, mas sao todas do corpo da soleira (nome do furo ja diz) -- as 3
    # amostras de referencia nao tem associacao espacial/litologica, ficam
    # de fora da contagem por corpo.
    for g in registros_geoq:
        litologia = g["litologia"]
        if not litologia and g["origem"] and g["origem"].startswith("Furo"):
            litologia = "sill_diabasio"
        if litologia not in resumo or not g["classificacao_ti"]:
            continue
        if g["classificacao_ti"] == "Alto-Ti":
            resumo[litologia]["alto"] += 1
        elif g["classificacao_ti"] == "Baixo-Ti":
            resumo[litologia]["baixo"] += 1
    return list(resumo.values())


# ======================================================================
# 1d. apanhado estatistico dos elementos quimicos (oxidos maiores + tracos)
#     -- pro quadro embaixo da coluna de graficos.
# ======================================================================
ELEMENTOS_OXIDOS = [
    ("SiO2", "SiO₂", "%"), ("TiO2", "TiO₂", "%"), ("Al2O3", "Al₂O₃", "%"), ("Fe2O3T", "Fe₂O₃(t)", "%"),
    ("MgO", "MgO", "%"), ("CaO", "CaO", "%"), ("Na2O", "Na₂O", "%"), ("K2O", "K₂O", "%"),
    ("MnO", "MnO", "%"), ("P2O5", "P₂O₅", "%"),
]
ELEMENTOS_TRACOS = [
    ("Sr", "Sr", "ppm"), ("Zr", "Zr", "ppm"), ("Y", "Y", "ppm"),
    ("Nb", "Nb", "ppm"), ("Rb", "Rb", "ppm"), ("Ba", "Ba", "ppm"),
]


def calcular_estatisticas_elementos():
    if not GEOQUIMICA_CSV.exists():
        return [], []
    df = pd.read_csv(GEOQUIMICA_CSV)

    tem_classificacao = "classificacao_ti" in df.columns

    def montar(lista_elementos):
        stats = []
        for coluna, label, unidade in lista_elementos:
            if coluna not in df.columns:
                continue
            serie = df[coluna].dropna()
            if serie.empty:
                continue
            media_alto = media_baixo = None
            if tem_classificacao:
                serie_alto = df.loc[df["classificacao_ti"] == "Alto-Ti", coluna].dropna()
                serie_baixo = df.loc[df["classificacao_ti"] == "Baixo-Ti", coluna].dropna()
                if not serie_alto.empty:
                    media_alto = float(serie_alto.mean())
                if not serie_baixo.empty:
                    media_baixo = float(serie_baixo.mean())
            stats.append({
                "label": label, "unidade": unidade, "n": int(serie.count()),
                "min": float(serie.min()), "media": float(serie.mean()),
                "mediana": float(serie.median()), "max": float(serie.max()),
                "media_alto": media_alto, "media_baixo": media_baixo,
            })
        return stats

    return montar(ELEMENTOS_OXIDOS), montar(ELEMENTOS_TRACOS)


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
NOTA_COM_DADO = "41 amostras (QMC_TAIO_TODOS) · campos de referência digitalizados aproximados (Fontoura, TCC 2024 · Peate et al. 1997)"


def trace_pontos_geoq(pontos, x_key, y_key, nome="Amostras Taió"):
    """Scatter dos pontos reais de geoquimica coloridos por Alto/Baixo-Ti,
    com hover mostrando a amostra e o ponto de campo associado (quando ha).
    customdata carrega o id (GEOQ-<amostra>) pra clique sincronizar com a
    lista/mapa via selecionarPorId no JS."""
    xs, ys, cores, textos, ids = [], [], [], [], []
    for p in pontos:
        x, y = p.get(x_key), p.get(y_key)
        if x is None or y is None:
            continue
        xs.append(x)
        ys.append(y)
        cores.append(p["cor"])
        rotulo = p["amostra"] + (f" ({NOME_PONTO.get(p['ponto_id'], p['ponto_id'])})" if p.get("ponto_id") and NOME_PONTO.get(p['ponto_id'], p['ponto_id']) != p["amostra"] else "")
        textos.append(rotulo)
        ids.append(p["id_display"])
    return go.Scatter(
        x=xs, y=ys, mode="markers", name=nome, text=textos, hoverinfo="text", customdata=ids,
        marker=dict(size=8, color=cores, line=dict(width=1, color="#1B1F2E")),
        showlegend=False,
    )


# ======================================================================
# 3. diagramas geoquimicos padrao (TAS, AFM, Shand's Index) -- limites bem
#    conhecidos da literatura (Le Bas et al. 1986 / Irvine & Baragar 1971),
#    simplificados (linhas principais, sem todas as subdivisoes finas) ja
#    que ainda nao ha ponto real do Taio pra classificar contra eles.
# ======================================================================
def montar_tas(pontos=None):
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
    if pontos:
        fig.add_trace(trace_pontos_geoq(pontos, "tas_x", "tas_y"))
    fig.update_xaxes(title_text="SiO₂", range=[38, 80])
    fig.update_yaxes(title_text="Na₂O + K₂O", range=[0, 16])
    return tema_grafico(fig, "TAS — Total Álcalis vs Sílica", altura=380, nota=NOTA_COM_DADO if pontos else NOTA_SEM_DADO)


def montar_afm(pontos=None):
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
    if pontos:
        xs_a, xs_b, xs_c, cores, textos, ids = [], [], [], [], [], []
        for p in pontos:
            if p.get("afm_a_F") is None:
                continue
            xs_a.append(p["afm_a_F"]); xs_b.append(p["afm_b_A"]); xs_c.append(p["afm_c_M"])
            cores.append(p["cor"])
            textos.append(p["amostra"] + (f" ({NOME_PONTO.get(p['ponto_id'], p['ponto_id'])})" if p.get("ponto_id") and NOME_PONTO.get(p['ponto_id'], p['ponto_id']) != p["amostra"] else ""))
            ids.append(p["id_display"])
        fig.add_trace(go.Scatterternary(
            a=xs_a, b=xs_b, c=xs_c, mode="markers", text=textos, hoverinfo="text", customdata=ids,
            name="Amostras Taió", marker=dict(size=8, color=cores, line=dict(width=1, color="#1B1F2E")), showlegend=False,
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
    return tema_grafico(fig, "AFM — Álcalis / FeOt / MgO", altura=380, nota=NOTA_COM_DADO if pontos else NOTA_SEM_DADO)


def montar_shand(pontos=None):
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
    if pontos:
        fig.add_trace(trace_pontos_geoq(pontos, "shand_acnk", "shand_ank"))
    fig.update_xaxes(title_text="A/CNK", range=[0, 1.6])
    fig.update_yaxes(title_text="A/NK", range=[0, 7])
    return tema_grafico(fig, "Shand's Index — A/NK × A/CNK", nota=NOTA_COM_DADO if pontos else NOTA_SEM_DADO)


def registro_lista_geoq(g):
    """Converte uma amostra de geoquimica bruta pro mesmo formato dos
    registros de campo, pra poder entrar na lista/DADOS_POR_ID e ficar
    clicavel igual (selecionarPorId funciona pra qualquer id presente ali,
    nao importa a origem)."""
    r = {
        "id": g["id"],
        "nome": g["amostra"],
        "litologia": g["litologia"] or "—",
        "classificacao_ti": g["classificacao_ti"],
        "tipo_ponto": f"Amostra ({g['origem']})",
        "qualidade": "",
        "x": g["utm_e"], "y": g["utm_n"],
        "descricao": f"Prof.: {g['profundidade_m']:.1f} m" if g["profundidade_m"] is not None else "",
        "tem_geoquimica": True,
        "cor_mapa": g["cor"],
        "geoq": {k: g[k] for k in ("SiO2", "TiO2", "Al2O3", "FeOT", "CaO", "MgO", "K2O", "Na2O", "P2O5",
                                    "Sr", "Zr", "Y", "Ti_Zr", "Zr_Y", "amostra")},
    }
    if g["lat"] is not None and g["lon"] is not None:
        r["lat"], r["lon"] = g["lat"], g["lon"]
    return r


# ======================================================================
# Critérios de Peate et al. (1992), matrizes de conformidade, ETR (manto primitivo),
# discriminação tectônica (Ti/V) e Harker -- todos a partir da planilha QMC
# (QMC Taio + Furos.xlsx), classificação Alto/Baixo-Ti vinda do CSV do dashboard.
# Janelas dos 12 critérios copiadas de Peate et al. (1992) -- NAO alterar.
# ======================================================================
QMC_XLSX = GEOQUIMICA_CSV.parent / "QMC Taio + Furos.xlsx"

CRITERIOS_PEATE = {
    'Urubici':      {'SiO2': (49, None), 'TiO2': (3.3, None), 'P2O5': (0.45, None), 'Fe2O3t': (None, 14.5),
                     'Sr': (550, None), 'Ba': (500, None), 'Zr': (250, None),
                     'TiZr': (57, None), 'TiY': (500, None), 'ZrY': (6.5, None), 'SrY': (14, None), 'BaY': (14, None)},
    'Pitanga':      {'SiO2': (47, None), 'TiO2': (2.8, None), 'P2O5': (0.35, None), 'Fe2O3t': (12.5, 18),
                     'Sr': (350, None), 'Ba': (200, None), 'Zr': (200, None),
                     'TiZr': (60, None), 'TiY': (350, None), 'ZrY': (5.5, None), 'SrY': (8, None), 'BaY': (9, None)},
    'Paranapanema': {'SiO2': (48, 53), 'TiO2': (1.7, 3.2), 'P2O5': (0.2, 0.8), 'Fe2O3t': (12.5, 17),
                     'Sr': (200, 450), 'Ba': (200, 650), 'Zr': (120, 250),
                     'TiZr': (65, None), 'TiY': (350, None), 'ZrY': (4.0, 7.0), 'SrY': (4.5, 15), 'BaY': (5, 19)},
    'Ribeira':      {'SiO2': (49, 52), 'TiO2': (1.5, 2.3), 'P2O5': (0.15, 0.50), 'Fe2O3t': (12, 16),
                     'Sr': (200, 375), 'Ba': (200, 600), 'Zr': (100, 200),
                     'TiZr': (65, None), 'TiY': (300, None), 'ZrY': (3.5, 7.0), 'SrY': (5, 17), 'BaY': (6, 19)},
    'Esmeralda':    {'SiO2': (48, 55), 'TiO2': (1.1, 2.3), 'P2O5': (0.1, 0.35), 'Fe2O3t': (12, 17),
                     'Sr': (None, 250), 'Ba': (90, 400), 'Zr': (65, 210),
                     'TiZr': (60, None), 'TiY': (None, 330), 'ZrY': (2.0, 5.0), 'SrY': (None, 9), 'BaY': (None, 12)},
    'Gramado':      {'SiO2': (49, 60), 'TiO2': (0.7, 2.0), 'P2O5': (0.05, 0.40), 'Fe2O3t': (9, 16),
                     'Sr': (140, 400), 'Ba': (100, 700), 'Zr': (65, 275),
                     'TiZr': (None, 70), 'TiY': (None, 330), 'ZrY': (3.5, 6.5), 'SrY': (None, 13), 'BaY': (None, 19)},
}
ORDEM_TIPOS = ['Urubici', 'Pitanga', 'Paranapanema', 'Ribeira', 'Esmeralda', 'Gramado']
ORDEM_CRITERIOS = ['SiO2', 'TiO2', 'P2O5', 'Fe2O3t', 'Sr', 'Ba', 'Zr', 'TiZr', 'TiY', 'ZrY', 'SrY', 'BaY']
TIPOS_GRUPO = {"Alto-Ti": ['Urubici', 'Pitanga', 'Paranapanema'], "Baixo-Ti": ['Ribeira', 'Esmeralda', 'Gramado']}
# tipos de alto-Ti em laranjas, de baixo-Ti em azuis (paleta do artigo)
COR_TIPO = {'Urubici': '#B35900', 'Pitanga': '#E8820C', 'Paranapanema': '#F6B76B',
            'Ribeira': '#2F6FA8', 'Esmeralda': '#5BA3D9', 'Gramado': '#1B4F7A'}
ROTULO_CRITERIO = {
    'SiO2': "SiO₂ (% peso)", 'TiO2': "TiO₂ (% peso)", 'P2O5': "P₂O₅ (% peso)", 'Fe2O3t': "Fe₂O₃(t) (% peso)",
    'Sr': "Sr (ppm)", 'Ba': "Ba (ppm)", 'Zr': "Zr (ppm)", 'TiZr': "Ti/Zr", 'TiY': "Ti/Y", 'ZrY': "Zr/Y",
    'SrY': "Sr/Y", 'BaY': "Ba/Y",
}
# manto primitivo (McDonough & Sun, 1995), ppm
ETR_ELEMENTOS = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
MANTO_PRIMITIVO = {'La': 0.648, 'Ce': 1.675, 'Pr': 0.254, 'Nd': 1.250, 'Sm': 0.406, 'Eu': 0.154, 'Gd': 0.544,
                   'Tb': 0.099, 'Dy': 0.674, 'Ho': 0.149, 'Er': 0.438, 'Tm': 0.068, 'Yb': 0.441, 'Lu': 0.0675}
HARKER_Y = [('TiO2', "TiO₂"), ('Al2O3', "Al₂O₃"), ('Fe2O3t', "Fe₂O₃(t)"), ('MgO', "MgO"), ('CaO', "CaO"),
            ('Na2O', "Na₂O"), ('K2O', "K₂O"), ('P2O5', "P₂O₅")]
NOTA_PEATE = "Peate et al. (1992) · faixas = janelas de cada magma-tipo · QMC_TAIO_TODOS"


def carregar_dados_peate(registros_geoq):
    """Planilha QMC + campos derivados (formulas da metodologia) + Alto/Baixo-Ti e
    id do dashboard (pra clique sincronizar com lista/mapa)."""
    df = pd.read_excel(QMC_XLSX, sheet_name="Planilha1")
    df = df[df["SAMPLE NAME"].notna()].copy()
    df["amostra"] = df["SAMPLE NAME"].astype(str).str.strip()
    por_nome = {r["amostra"]: r for r in registros_geoq}
    df = df[df["amostra"].isin(por_nome)].copy()
    df["classe"] = df["amostra"].map(lambda a: por_nome[a]["classificacao_ti"])
    df["id"] = df["amostra"].map(lambda a: por_nome[a]["id_display"])
    df["cor"] = df["amostra"].map(lambda a: por_nome[a]["cor"])
    df["SiO2"], df["TiO2"], df["Al2O3"] = df["SiO2_pct"], df["TiO2_pct"], df["Al2O3"]
    df["MgO"], df["CaO"], df["Na2O"], df["K2O"], df["P2O5"] = df["MgO_pct"], df["CaO_pct"], df["Na2O_pct"], df["K2O_pct"], df["P2O5_pct"]
    df["Fe2O3t"] = df["FeOT_pct"] * 1.1113
    df["Ti_ppm"] = df["TiO2_pct"] * 5990.0
    df["Sr"], df["Ba"], df["Zr"], df["Y"] = df["Sr_ppm"], df["Ba_ppm"], df["Zr_ppm"], df["Y_ppm"]
    df["TiZr"], df["TiY"] = df["Ti_ppm"] / df["Zr"], df["Ti_ppm"] / df["Y"]
    df["ZrY"], df["SrY"], df["BaY"] = df["Zr"] / df["Y"], df["Sr"] / df["Y"], df["Ba"] / df["Y"]
    return df.reset_index(drop=True)


def _dentro(valor, janela):
    lo, hi = janela
    return (lo is None or valor >= lo) and (hi is None or valor <= hi)


def pontuacao_peate(df):
    """Pontos 0-12 por amostra x tipo: 1 por criterio dentro da janela do tipo."""
    return pd.DataFrame({
        tipo: [sum(_dentro(row[c], CRITERIOS_PEATE[tipo][c]) for c in ORDEM_CRITERIOS) for _, row in df.iterrows()]
        for tipo in ORDEM_TIPOS
    }, index=df["amostra"])


def trace_amostras(df, x, y, **extra):
    return go.Scatter(
        x=x, y=y, mode="markers", name="Amostras Taió", text=df["amostra"], customdata=df["id"], hoverinfo="text+y",
        marker=dict(size=8, color=df["cor"], line=dict(width=1, color="#071729")), showlegend=False, **extra,
    )


def montar_criterio_peate(df, crit):
    """Um criterio: valor de cada amostra (ordenadas por grupo e valor) com as janelas dos 6 tipos ao fundo."""
    d = df.sort_values(["classe", crit]).reset_index(drop=True)
    vals = d[crit].astype(float)
    finitos = [v for t in ORDEM_TIPOS for v in CRITERIOS_PEATE[t][crit] if v is not None]
    lo = min([vals.min()] + finitos)
    hi = max([vals.max()] + finitos)
    folga = (hi - lo) * 0.06 or 1.0
    y0, y1 = lo - folga, hi + folga
    fig = go.Figure()
    shapes = []
    for tipo in ORDEM_TIPOS:
        jmin, jmax = CRITERIOS_PEATE[tipo][crit]
        a, b = (y0 if jmin is None else jmin), (y1 if jmax is None else jmax)
        cor = COR_TIPO[tipo]
        shapes.append(dict(type="rect", xref="paper", x0=0, x1=1, yref="y", y0=a, y1=b, fillcolor=cor, opacity=0.14,
                           line=dict(width=0), layer="below"))
        for yy in (jmin, jmax):
            if yy is not None:
                shapes.append(dict(type="line", xref="paper", x0=0, x1=1, yref="y", y0=yy, y1=yy,
                                   line=dict(color=cor, width=1, dash="dot"), layer="below"))
        # trace invisivel so pra hover da janela
        fig.add_trace(go.Scatter(x=[0.5], y=[(a + b) / 2], mode="markers", marker=dict(size=0.1, color=cor), xaxis="x2",
                                 showlegend=False, hoverinfo="text", text=[f"{tipo}: {jmin if jmin is not None else '—'} a {jmax if jmax is not None else '—'}"]))
    fig.add_trace(trace_amostras(d, list(range(len(d))), vals))
    fig.update_layout(shapes=shapes, xaxis2=dict(overlaying="x", range=[0, 1], visible=False))
    fig.update_xaxes(showticklabels=False, range=[-1, len(d)], title_text="amostras (Alto-Ti → Baixo-Ti)")
    fig.update_yaxes(title_text=ROTULO_CRITERIO[crit], range=[y0, y1])
    tema_grafico(fig, ROTULO_CRITERIO[crit], altura=250)
    fig.update_layout(margin=dict(l=48, r=10, t=36, b=42), title=dict(font=dict(size=12)), hovermode="closest")
    return fig


def montar_matriz_conformidade(df, pont, grupo):
    """Heatmap amostras x tipos do proprio grupo (nao reclassifica Alto-Ti contra tipos de Baixo-Ti)."""
    tipos = TIPOS_GRUPO[grupo]
    sel = df[df["classe"] == grupo]
    sub = pont.loc[sel["amostra"], tipos]
    ordem = sub.assign(_m=sub.max(axis=1), _s=sub.sum(axis=1)).sort_values(["_m", "_s"], ascending=False).index
    sub = sub.loc[ordem]
    ids = sel.set_index("amostra").loc[ordem, "id"].tolist()
    cor_max = "#e8820c" if grupo == "Alto-Ti" else "#2f6fa8"
    fig = go.Figure(go.Heatmap(
        z=sub.values, x=tipos, y=list(sub.index), customdata=[[i] * len(tipos) for i in ids],
        colorscale=[[0, "#ffffff"], [1, cor_max]], zmin=0, zmax=12, xgap=2, ygap=2,
        text=sub.values, texttemplate="%{text}", textfont=dict(size=10, color="#111"),
        colorbar=dict(title="/12", thickness=10, len=0.6), hovertemplate="%{y} · %{x}: %{z}/12<extra></extra>",
    ))
    fig.update_yaxes(autorange="reversed", type="category", tickfont=dict(size=9))
    fig.update_xaxes(side="top", type="category")
    tema_grafico(fig, f"Conformidade com os 6 critérios de Peate — {grupo} ({len(sub)} amostras)", altura=max(260, 17 * len(sub) + 90),
                 nota="pontos = nº de critérios (0–12) dentro da janela de cada magma-tipo do grupo · ordenado pela melhor pontuação")
    fig.update_layout(margin=dict(l=62, r=20, t=70, b=70))
    return fig


def montar_etr(df):
    completos = df.dropna(subset=[f"{e}_ppm" for e in ETR_ELEMENTOS])
    fig = go.Figure()
    vistos = set()
    for _, r in completos.iterrows():
        y = [r[f"{e}_ppm"] / MANTO_PRIMITIVO[e] for e in ETR_ELEMENTOS]
        fig.add_trace(go.Scatter(
            x=ETR_ELEMENTOS, y=y, mode="lines+markers", name=r["classe"], legendgroup=r["classe"],
            showlegend=r["classe"] not in vistos, text=[r["amostra"]] * len(y), customdata=[r["id"]] * len(y),
            hovertemplate="%{text} · %{x}: %{y:.1f}×MP<extra></extra>", line=dict(color=r["cor"], width=1.3),
            marker=dict(size=4, color=r["cor"]), opacity=0.8,
        ))
        vistos.add(r["classe"])
    fig.update_yaxes(type="log", title_text="Amostra / manto primitivo")
    fig.update_xaxes(title_text="Elementos terras raras")
    tema_grafico(fig, f"Padrões de ETR normalizados ao manto primitivo ({len(completos)} amostras)", altura=360,
                 nota="normalização: manto primitivo, McDonough & Sun (1995) · só amostras com os 14 ETR analisados")
    return fig


def montar_ti_v(df):
    """Discriminação tectônica Ti/V (Shervais, 1982): Ti/1000 x V."""
    fig = go.Figure()
    vmax = max(650, float(df["V_ppm"].max()) * 1.1)
    for razao, rot in ((10, "10"), (20, "20"), (50, "50"), (100, "100")):
        fig.add_trace(go.Scatter(x=[0, vmax], y=[0, razao * vmax / 1000], mode="lines", showlegend=False, hoverinfo="skip",
                                 line=dict(color="#8A94A6", width=1, dash="dot")))
    ymax = max(df["Ti_ppm"].max() / 1000 * 1.15, 50 * vmax / 1000 * 0.35)
    campos = [(0.8 * vmax, 0.8 * vmax * 15 / 1000, "Ti/V < 20 · arco"), (0.8 * vmax, 0.8 * vmax * 35 / 1000, "20–50 · MORB / retroarco"),
              (0.8 * vmax, 0.8 * vmax * 75 / 1000, "> 50 · OIB / CFB")]
    ann = [dict(x=x, y=y, text=t, showarrow=False, font=dict(size=9, color="#8A94A6"), xanchor="right") for x, y, t in campos]
    for razao in (10, 20, 50, 100):
        ann.append(dict(x=vmax, y=razao * vmax / 1000, text=str(razao), showarrow=False, font=dict(size=9, color="#8A94A6"),
                        xanchor="left", yanchor="middle"))
    fig.add_trace(trace_amostras(df, df["V_ppm"], df["Ti_ppm"] / 1000))
    fig.update_layout(annotations=ann)
    fig.update_xaxes(title_text="V (ppm)", range=[0, vmax])
    fig.update_yaxes(title_text="Ti / 1000 (ppm)", range=[0, ymax])
    return tema_grafico(fig, "Discriminação tectônica — Ti/V (Shervais, 1982)", altura=360,
                        nota="linhas = razões Ti/V 10, 20, 50 e 100 · Ti = TiO₂ × 5990")


def montar_harker(df, col, rotulo):
    fig = go.Figure(trace_amostras(df, df["SiO2"], df[col]))
    fig.update_xaxes(title_text="SiO₂ (% peso)")
    fig.update_yaxes(title_text=f"{rotulo} (% peso)")
    tema_grafico(fig, f"Harker — {rotulo} × SiO₂", altura=250)
    fig.update_layout(margin=dict(l=48, r=10, t=36, b=42), title=dict(font=dict(size=12)))
    return fig


def montar_blocos_peate(registros_geoq):
    """Devolve (lista de (div_id, fig), dict de secoes) pro HTML do dashboard."""
    df = carregar_dados_peate(registros_geoq)
    pont = pontuacao_peate(df)
    figs = {}
    for c in ORDEM_CRITERIOS:
        figs[f"grafico-peate-{c.lower()}"] = montar_criterio_peate(df, c)
    figs["grafico-matriz-alto"] = montar_matriz_conformidade(df, pont, "Alto-Ti")
    figs["grafico-matriz-baixo"] = montar_matriz_conformidade(df, pont, "Baixo-Ti")
    figs["grafico-etr"] = montar_etr(df)
    figs["grafico-ti-v"] = montar_ti_v(df)
    for col, rot in HARKER_Y:
        figs[f"grafico-harker-{col.lower()}"] = montar_harker(df, col, rot)
    print(f"  Peate/ETR/Ti-V/Harker: {len(df)} amostras, {len(figs)} gráficos")
    return figs


def main():
    print("Carregando pontos de campo (Taió)...")
    registros_campo = carregar_campo()
    NOME_PONTO.update({r["ponto_id"]: r["nome"] for r in registros_campo})
    print(f"  {len(registros_campo)} pontos")

    print("Carregando geoquímica bruta (QMC_TAIO_TODOS)...")
    registros_geoq = carregar_geoquimica()
    print(f"  {len(registros_geoq)} amostras ({sum(1 for r in registros_geoq if r['lat'] is not None)} com coordenadas)")
    geoq_por_ponto = {r["ponto_id"]: r for r in registros_geoq if r["ponto_id"]}
    # unifica a classificacao Ti: a geoquimica bruta real (QMC_TAIO_TODOS) e a
    # fonte de verdade quando existe pro ponto -- sobrescreve o campo antigo
    # do catalogo (ti_geoquimico, pontos_unificados_completo.gpkg), que podia
    # estar desatualizado/faltando em relacao as 41 amostras reais. Sem isso
    # o grafico "Classificação de Ti" ficava quase todo "Sem dado" mesmo pros
    # pontos que ja tem amostra real cruzada.
    MAPA_TI_GEOQ = {"Alto-Ti": "Alto", "Baixo-Ti": "Baixo"}
    for r in registros_campo:
        g = geoq_por_ponto.get(r["ponto_id"])
        if g:
            r["geoq"] = {k: g[k] for k in ("SiO2", "TiO2", "Al2O3", "FeOT", "CaO", "MgO", "K2O", "Na2O", "P2O5",
                                            "Sr", "Zr", "Y", "Ti_Zr", "Zr_Y", "amostra")}
            if g["classificacao_ti"] in MAPA_TI_GEOQ:
                r["classificacao_ti"] = MAPA_TI_GEOQ[g["classificacao_ti"]]

    print("Carregando resumo dos corpos intrusivos (área + classificação)...")
    resumo_corpos = carregar_resumo_corpos(registros_geoq)
    for corpo in resumo_corpos:
        print(f"  {corpo['nome']}: {corpo['area_km2']:.1f} km² ({corpo['n_poligonos']} polígonos) — "
              f"{corpo['alto']} Alto-Ti, {corpo['baixo']} Baixo-Ti")

    print("Calculando estatísticas dos elementos químicos...")
    estat_oxidos, estat_tracos = calcular_estatisticas_elementos()
    print(f"  {len(estat_oxidos)} óxidos, {len(estat_tracos)} traços")

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

    # amostras de geoquimica bruta sem lat/lon proprio na planilha mas com
    # ponto_id (ex.: ITC-44B) usam a coordenada do ponto de campo
    # correspondente (ja convertida pra WGS84 acima) -- so as 3 "amostras de
    # referencia" (sem ponto/corpo associado, coordenada real desconhecida)
    # ficam mesmo de fora do mapa.
    campo_por_nome = {r["ponto_id"]: r for r in registros_campo}
    for r in registros_geoq:
        if r["lat"] is None and r["ponto_id"]:
            c = campo_por_nome.get(r["ponto_id"])
            if c and c.get("lat") is not None:
                r["lat"], r["lon"] = c["lat"], c["lon"]

    # amostras de geoquimica bruta com coordenada (campo + furo -- furo tem
    # as 15 amostras na mesma coordenada de boca de furo, so a profundidade
    # muda; ainda assim cada uma vira um marcador, o popup mostra a profundidade)
    geoq_coord = [r for r in registros_geoq if r["lat"] is not None]
    geojson_geoq = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature", "geometry": {"type": "Point", "coordinates": [r["lon"], r["lat"]]},
            "properties": {
                "id": r["id_display"],
                "cor": r["cor"],
                "popup": (
                    f"<b>{r['amostra']}</b> ({r['classificacao_ti']})<br>{r['origem']}"
                    + (f"<br>Prof.: {r['profundidade_m']:.1f} m" if r["profundidade_m"] is not None else "")
                    + f"<br>SiO₂ {r['SiO2']:.1f}% · TiO₂ {r['TiO2']:.2f}% · MgO {r['MgO']:.2f}% · P₂O₅ {r['P2O5']:.2f}%"
                ),
            },
        } for r in geoq_coord],
    }

    # mapa geologico atualizado (6 formacoes) -- mesma camada usada no
    # webmap dedicado (gerar_webmap_taio.py), pedido explicito do usuario
    # ("faltou o mapa geológico junto").
    CORES_LITOLOGIA_MAPA = dict(zip(NOMES_CAMADAS, CORES_CAMADAS))
    CORES_LITOLOGIA_MAPA["Depósito quaternário"] = COR_QUATERNARIO
    CORES_LITOLOGIA_MAPA["Soleira"] = COR_SILL
    CORES_LITOLOGIA_MAPA["Dique"] = COR_DIQUE
    geojson_formacoes = None
    if LITOLOGIA_ATUALIZADA.exists():
        gdf_formacoes = gpd.read_file(LITOLOGIA_ATUALIZADA).to_crs(4326)
        gdf_formacoes["cor"] = gdf_formacoes["formacao"].map(CORES_LITOLOGIA_MAPA).fillna("#CCCCCC")
        gdf_formacoes["popup"] = gdf_formacoes["formacao"]
        geojson_formacoes = json.loads(gdf_formacoes[["formacao", "cor", "popup", "geometry"]].to_json())
        print(f"  mapa geológico: {len(gdf_formacoes)} polígonos (incl. sill/dique)")

    # ------------------------------------------------------------
    # diagramas geoquimicos -- todos so com campos de referencia (sem
    # pontos, ver nota no topo do arquivo)
    # ------------------------------------------------------------
    fig_tas = montar_tas(pontos=registros_geoq)
    fig_afm = montar_afm(pontos=registros_geoq)
    fig_shand = montar_shand(pontos=registros_geoq)
    figs_peate = montar_blocos_peate(registros_geoq)

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
    # amostras de furo/referencia (QMC_TAIO_TODOS) nao correspondem a nenhum
    # ponto do catalogo de campo (sem ponto_id) -- somadas aqui direto, senao
    # ficariam de fora da contagem unificada de Alto/Baixo-Ti.
    for g in registros_geoq:
        if not g["ponto_id"]:
            contagem_ti_taio[MAPA_TI_GEOQ.get(g["classificacao_ti"], "Sem dado")] += 1
    ordem_ti = ["Alto", "Baixo", "Sem dado"]
    itens_ti = [(k, contagem_ti_taio.get(k, 0)) for k in ordem_ti if contagem_ti_taio.get(k, 0)]
    fig_ti_taio = tema_grafico(go.Figure(go.Bar(
        x=[k for k, _ in itens_ti], y=[v for _, v in itens_ti],
        marker=dict(color=["#E67E22", "#2E86C1", COR_TI_INDEFINIDO][:len(itens_ti)]),
    )), "Classificação de Ti")

    print("Montando HTML final...")
    # a biblioteca plotly.js so precisa ser embutida uma vez -- na PRIMEIRA
    # figura que aparece na pagina.
    html_tas = pio.to_html(fig_tas, full_html=False, include_plotlyjs=True, div_id="grafico-tas", config={"responsive": True})
    html_afm = pio.to_html(fig_afm, full_html=False, include_plotlyjs=False, div_id="grafico-afm", config={"responsive": True})
    html_shand = pio.to_html(fig_shand, full_html=False, include_plotlyjs=False, div_id="grafico-shand", config={"responsive": True})
    html_peate = {k: pio.to_html(f, full_html=False, include_plotlyjs=False, div_id=k, config={"responsive": True})
                  for k, f in figs_peate.items()}

    def cartoes(chaves):
        return "".join(f'<div class="grafico-card">{html_peate[k]}</div>' for k in chaves)

    legenda_tipos = "".join(f'<span><i style="background:{COR_TIPO[tp]}"></i>{tp}</span>' for tp in ORDEM_TIPOS)
    secao_peate = (
        '<h3 class="titulo-secao">Critérios de Peate et al. (1992) — 12 critérios</h3>'
        f'<div class="legenda-tipos">{legenda_tipos}<em>faixas = janela de cada magma-tipo · pontos = amostras</em></div>'
        f'<div class="grade-2">{cartoes([f"grafico-peate-{c.lower()}" for c in ORDEM_CRITERIOS])}</div>')
    secao_matrizes = (
        '<h3 class="titulo-secao">Matrizes de conformidade (Peate et al., 1992)</h3>'
        + cartoes(["grafico-matriz-alto", "grafico-matriz-baixo"]))
    secao_etr_tect = (
        '<h3 class="titulo-secao">ETR e discriminação tectônica</h3>' + cartoes(["grafico-etr", "grafico-ti-v"]))
    secao_harker = (
        '<h3 class="titulo-secao">Diagramas de Harker</h3>'
        f'<div class="grade-2">{cartoes([f"grafico-harker-{c.lower()}" for c, _ in HARKER_Y])}</div>')
    ids_novos_js = ", ".join(f"document.getElementById('{k}')" for k in figs_peate)
    html_lito_taio = pio.to_html(fig_lito_taio, full_html=False, include_plotlyjs=False, div_id="grafico-lito-taio", config={"responsive": True})
    html_ti_taio = pio.to_html(fig_ti_taio, full_html=False, include_plotlyjs=False, div_id="grafico-ti-taio", config={"responsive": True})

    # so amostras SEM ponto de campo correspondente (furo/referencia) viram
    # linha propria na lista -- as que ja tem ponto_id foram mescladas dentro
    # do registro de campo (r["geoq"]) la em cima, mostrar as duas seria
    # duplicar a mesma amostra fisica com dois nomes diferentes na lista.
    registros_geoq_sem_campo = [g for g in registros_geoq if not g["ponto_id"]]
    registros_geoq_lista = [registro_lista_geoq(g) for g in registros_geoq_sem_campo]
    dados_js = json.dumps(registros_campo + registros_geoq_lista, ensure_ascii=False)
    print("Calculando colunas estratigráficas estimadas...")
    colunas, max_total_coluna = calcular_colunas(registros_campo + registros_geoq_lista)
    print(f"  {len(colunas)} colunas (profundidade máxima {max_total_coluna:.0f} m)")
    colunas_js = json.dumps(colunas, ensure_ascii=False)
    logo_b64 = logo_base64()

    linhas_tabela = []
    for r in registros_campo:
        linhas_tabela.append(f"""
        <tr class="linha-dado" data-id="{r['id']}"
            data-busca="{(r['nome'] + ' ' + r['ponto_id'] + ' ' + (r['litologia'] or '')).lower()}">
            <td><span class="dot-cor" style="background:{r['cor_mapa']}"></span><span class="nome-ponto">{r['nome']}</span></td>
            <td>{r['litologia'] or '—'}</td>
            <td>{r['classificacao_ti'] or '—'}</td>
        </tr>""")
    for g in registros_geoq_sem_campo:
        busca = f"{g['amostra']} {g['origem']} amostra geoquimica {g['litologia'] or ''} {g['ponto_id'] or ''}".lower()
        linhas_tabela.append(f"""
        <tr class="linha-dado linha-geoq" data-id="{g['id']}" data-busca="{busca}">
            <td><span class="dot-cor" style="background:{g['cor']}"></span><span class="nome-ponto">{g['amostra']}</span></td>
            <td>{g['litologia'] or g['origem']}</td>
            <td>{g['classificacao_ti'] or '—'}</td>
        </tr>""")
    tabela_html = "".join(linhas_tabela)

    linhas_resumo = []
    for corpo in resumo_corpos:
        total_ti = corpo["alto"] + corpo["baixo"]
        pct_alto = (corpo["alto"] / total_ti * 100) if total_ti else 0
        pct_baixo = 100 - pct_alto if total_ti else 0
        barra = (
            f'<div class="barra-ti"><div class="seg-alto" style="width:{pct_alto:.1f}%"></div>'
            f'<div class="seg-baixo" style="width:{pct_baixo:.1f}%"></div></div>'
            f'<div class="legenda-barra">{corpo["alto"]} Alto · {corpo["baixo"]} Baixo</div>'
            if total_ti else '<div class="legenda-barra">Sem amostra classificada</div>'
        )
        volume_txt = f"{corpo['volume_km3']:.2f} km³" if corpo["volume_km3"] is not None else "—"
        linhas_resumo.append(f"""
        <tr>
            <td>{corpo['nome']}</td>
            <td>{corpo['area_km2']:.1f} km²</td>
            <td>{volume_txt}</td>
            <td>{corpo['n_poligonos']}</td>
            <td>{barra}</td>
        </tr>""")
    resumo_corpos_html = "".join(linhas_resumo)

    def montar_linhas_estatisticas(lista_stats):
        linhas = []
        for e in lista_stats:
            faixa = e["max"] - e["min"]
            pct_media = ((e["media"] - e["min"]) / faixa * 100) if faixa else 50
            pct_mediana = ((e["mediana"] - e["min"]) / faixa * 100) if faixa else 50
            casas = 2 if e["unidade"] == "%" else 1
            txt_alto = f"{e['media_alto']:.{casas}f}" if e["media_alto"] is not None else "—"
            txt_baixo = f"{e['media_baixo']:.{casas}f}" if e["media_baixo"] is not None else "—"
            linhas.append(f"""
            <tr>
                <td>{e['label']}</td>
                <td>{e['min']:.{casas}f}</td>
                <td>
                    <div class="faixa-elemento">
                        <div class="marca-mediana" style="left:{pct_mediana:.1f}%" title="Mediana: {e['mediana']:.{casas}f}"></div>
                        <div class="marca-media" style="left:{pct_media:.1f}%" title="Média: {e['media']:.{casas}f}"></div>
                    </div>
                </td>
                <td>{e['media']:.{casas}f} {e['unidade']}</td>
                <td style="color:{COR_TI_ALTO}">{txt_alto}</td>
                <td style="color:{COR_TI_BAIXO}">{txt_baixo}</td>
                <td>{e['max']:.{casas}f}</td>
            </tr>""")
        return "".join(linhas)

    estatisticas_oxidos_html = montar_linhas_estatisticas(estat_oxidos)
    estatisticas_tracos_html = montar_linhas_estatisticas(estat_tracos)

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
  .col-esquerda {{ display: flex; flex-direction: column; gap: 12px; min-height: 0; height: 100%; }}
  .painel {{ background: var(--painel); border: 1px solid {MARCA_ROXO}; border-radius: 8px; overflow: hidden; display: flex; flex-direction: column; }}
  .col-esquerda > .painel:first-child {{ flex: 1; min-height: 0; }}
  .painel-estatico {{ flex: 0 0 auto; }}
  .tabela-resumo {{ width: 100%; border-collapse: collapse; font-size: 11px; margin: 8px 14px; width: calc(100% - 28px); }}
  .tabela-resumo th {{ text-align: left; padding: 4px 6px; border-bottom: 1px solid {MARCA_ROXO}; opacity: 0.75; font-weight: 500; }}
  .tabela-resumo td {{ padding: 5px 6px; border-bottom: 1px solid var(--borda-fraca); vertical-align: middle; }}
  .tabela-resumo tr:last-child td {{ border-bottom: none; }}
  .barra-ti {{ display: flex; width: 100%; height: 8px; border-radius: 4px; overflow: hidden; background: var(--borda-fraca); margin-bottom: 3px; }}
  .barra-ti .seg-alto {{ background: {COR_TI_ALTO}; }}
  .barra-ti .seg-baixo {{ background: {COR_TI_BAIXO}; }}
  .legenda-barra {{ font-size: 10px; opacity: 0.65; }}
  .nota-resumo {{ font-size: 10px; opacity: 0.5; margin: 0 14px 10px 14px; line-height: 1.4; }}
  .subtitulo-estatistica {{ font-size: 11px; font-weight: 600; opacity: 0.75; margin: 10px 14px 0 14px; text-transform: uppercase; letter-spacing: 0.04em; }}
  .tabela-estatistica td:first-child {{ white-space: nowrap; }}
  .faixa-elemento {{ position: relative; width: 100%; min-width: 70px; height: 10px; }}
  .faixa-elemento::before {{ content: ''; position: absolute; top: 4px; left: 0; right: 0; height: 2px; background: var(--borda-fraca); border-radius: 1px; }}
  .marca-media {{ position: absolute; top: 1px; width: 8px; height: 8px; border-radius: 50%; background: {MARCA_ROXO}; transform: translateX(-50%); }}
  .marca-mediana {{ position: absolute; top: -1px; width: 2px; height: 12px; background: var(--texto); opacity: 0.55; transform: translateX(-50%); }}
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
  .dot-cor {{ display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 5px; }}
  .linha-geoq {{ border-left: 3px solid {MARCA_ROXO}; }}
  footer {{ text-align: center; padding: 8px; opacity: 0.55; font-size: 11px; }}
  #mapa-leaflet {{ flex: 1; }}
  .col-graficos {{ overflow-y: auto; padding: 10px; display: flex; flex-direction: column; gap: 12px; }}
  .titulo-secao {{ margin: 8px 2px 0; font-size: 13px; letter-spacing: .06em; text-transform: uppercase; color: {MARCA_ROXO}; flex-shrink: 0; }}
  .legenda-tipos {{ display: flex; flex-wrap: wrap; gap: 4px 12px; font-size: 11px; align-items: center; flex-shrink: 0; }}
  .legenda-tipos i {{ display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 4px; }}
  .legenda-tipos em {{ opacity: .6; }}
  .grade-2 {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 10px; flex-shrink: 0; }}
  .grade-2 .grafico-card {{ min-width: 0; overflow: hidden; }}
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
    background: {MARCA_ROXO_ESCURO}; color: {MARCA_CINZA_CLARO}; border: 1px solid {MARCA_ROXO}; border-radius: 8px;
    padding: 10px 12px; font-size: 12px; box-shadow: 0 4px 16px rgba(0,0,0,0.5); pointer-events: none;
  }}
  .coluna-corpo {{ display: flex; gap: 6px; padding: 12px 14px 14px 14px; align-items: flex-start; }}
  .coluna-eixo {{ position: relative; width: 38px; flex: none; font-size: 10px; opacity: 0.7; }}
  .coluna-eixo span {{ position: absolute; right: 4px; transform: translateY(-50%); }}
  .coluna-pilha {{ position: relative; width: 210px; flex: none; border: 1px solid var(--borda-fraca); }}
  .coluna-seg {{
    position: absolute; left: 0; right: 0; box-sizing: border-box; overflow: hidden;
    display: flex; flex-direction: column; align-items: center; justify-content: center;
    text-align: center; font-size: 11px; line-height: 1.25; border-bottom: 1px solid rgba(255,255,255,0.55);
  }}
  .coluna-lado {{ position: relative; width: 14px; flex: none; }}
  .coluna-dique {{ position: absolute; left: 3px; width: 8px; background: {COR_DIQUE}; border: 1px solid #000; }}
  .coluna-info {{ flex: 1; font-size: 11px; line-height: 1.5; opacity: 0.9; min-width: 0; }}
  .coluna-info .coluna-titulo {{ font-size: 13px; font-weight: 700; margin-bottom: 4px; }}
  .coluna-vazia {{ padding: 14px; font-size: 12px; opacity: 0.7; }}
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
  <div class="col-esquerda">
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
    <div class="painel painel-estatico">
      <h2>Corpos intrusivos — área, volume &amp; classificação</h2>
      <table class="tabela-resumo">
        <thead><tr><th>Corpo</th><th>Área</th><th>Volume</th><th>Polígonos</th><th>Classificação Ti</th></tr></thead>
        <tbody>{resumo_corpos_html}</tbody>
      </table>
      <p class="nota-resumo">Área digitalizada em planta (poligon_intrusiva.shp) · volume estimado (sill: topografia real menos topo da Serra Alta, piso de 27m · dique: espessura mediana de campo de 6m × área) · classificação Ti das amostras de geoquímica (furo contam como soleira; 3 amostras de referência sem corpo associado ficam de fora)</p>
    </div>
  </div>
  <div class="painel col-graficos">
    <div class="grafico-card painel-estatico">
      <h2>Coluna estratigráfica estimada (furo virtual)</h2>
      <div id="coluna-estrat"><div class="coluna-vazia">Selecione um ponto na lista ou no mapa para ver o empilhamento esperado ali.</div></div>
    </div>
    <div class="grafico-card">{html_tas}</div>
    <div class="grafico-card">{html_afm}</div>
    <div class="grafico-card">{html_shand}</div>
    {secao_peate}
    {secao_matrizes}
    {secao_etr_tect}
    {secao_harker}
    <div class="grafico-card">{html_lito_taio}</div>
    <div class="grafico-card">{html_ti_taio}</div>
    <div class="grafico-card painel-estatico">
      <h2>Apanhado estatístico — elementos químicos (41 amostras)</h2>
      <p class="subtitulo-estatistica">Óxidos maiores (% em peso)</p>
      <table class="tabela-resumo tabela-estatistica">
        <thead><tr><th>Elemento</th><th>Mín</th><th>Distribuição</th><th>Média</th><th style="color:{COR_TI_ALTO}">Média Alto-Ti</th><th style="color:{COR_TI_BAIXO}">Média Baixo-Ti</th><th>Máx</th></tr></thead>
        <tbody>{estatisticas_oxidos_html}</tbody>
      </table>
      <p class="subtitulo-estatistica">Elementos traço (ppm)</p>
      <table class="tabela-resumo tabela-estatistica">
        <thead><tr><th>Elemento</th><th>Mín</th><th>Distribuição</th><th>Média</th><th style="color:{COR_TI_ALTO}">Média Alto-Ti</th><th style="color:{COR_TI_BAIXO}">Média Baixo-Ti</th><th>Máx</th></tr></thead>
        <tbody>{estatisticas_tracos_html}</tbody>
      </table>
      <p class="nota-resumo">Faixa vai do mínimo ao máximo observado · <span style="color:{MARCA_ROXO}">●</span> média · <span style="opacity:0.6">▏</span> mediana · médias por classificação Ti calculadas só sobre as amostras de cada grupo</p>
    </div>
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
        document.getElementById('grafico-shand'),
        {ids_novos_js},
        document.getElementById('grafico-lito-taio'), document.getElementById('grafico-ti-taio'),
    ];

    // clique num ponto de amostra em qualquer diagrama sincroniza com a
    // lista/mapa (mesmo id GEOQ-<amostra> usado no customdata da trace).
    TODOS_GD.forEach(function(gd) {{
        if (!gd) return;
        gd.on('plotly_click', function(ev) {{
            var pt = ev.points && ev.points[0];
            if (!pt || pt.customdata === undefined || pt.customdata === null) return;
            selecionarPorId(pt.customdata, 'grafico');
        }});
    }});

    // indice (grafico, indice-da-trace, lista-de-ids) de todo grafico que
    // tem a trace de amostras reais (mesmo "name" usado em todos eles) --
    // usado pra destacar (aumentar/contornar) o ponto selecionado, e tirar
    // o destaque de qualquer outro grafico ao trocar de selecao.
    var GRAFICOS_PONTOS = [];
    TODOS_GD.forEach(function(gd) {{
        if (!gd || !gd.data) return;
        var idx = gd.data.findIndex(function(tr) {{ return tr.name === 'Amostras Taió'; }});
        if (idx === -1) return;
        GRAFICOS_PONTOS.push({{ gd: gd, idx: idx, customdata: (gd.data[idx].customdata || []).slice() }});
    }});

    function atualizarDestaqueGraficos(idSelecionado) {{
        GRAFICOS_PONTOS.forEach(function(g) {{
            var n = g.customdata.length;
            var tamanhos = new Array(n).fill(8);
            var largurasLinha = new Array(n).fill(1);
            var coresLinha = new Array(n).fill('#1B1F2E');
            var pos = g.customdata.indexOf(idSelecionado);
            if (pos !== -1) {{
                tamanhos[pos] = 14;
                largurasLinha[pos] = 3;
                coresLinha[pos] = '{MARCA_ROXO}';
            }}
            Plotly.restyle(g.gd, {{
                'marker.size': [tamanhos], 'marker.line.width': [largurasLinha], 'marker.line.color': [coresLinha],
            }}, [g.idx]);
        }});
    }}

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
    satelite.addTo(mapa);
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
    }});
    var geoqLayer = L.geoJSON({json.dumps(geojson_geoq, ensure_ascii=False)}, {{
        pointToLayer: function(f, latlng) {{
            return L.circleMarker(latlng, {{
                radius: 7, fillColor: f.properties.cor, color: '#1B1F2E', weight: 1.5,
                fillOpacity: 0.95, pane: 'markerPane',
            }});
        }},
        onEachFeature: function(f, layer) {{
            layer.bindPopup(f.properties.popup);
            layer.on('click', function() {{ selecionarPorId(f.properties.id, 'mapa'); }});
        }},
    }}).addTo(mapa);
    var overlaysMapa = {{ "Campo (Taió)": campoLayer, "Geoquímica (amostras brutas)": geoqLayer }};"""
    if geojson_formacoes is not None:
        html_final += f"""
    var formacoesLayer = L.geoJSON({json.dumps(geojson_formacoes, ensure_ascii=False)}, {{
        style: function(f) {{ return {{ color: '#000', weight: 0.5, fillColor: f.properties.cor, fillOpacity: 0.5 }}; }},
        onEachFeature: function(f, layer) {{ layer.bindPopup(f.properties.popup); }},
    }});
    overlaysMapa["Mapa geológico atualizado"] = formacoesLayer;"""
    html_final += f"""
    L.control.layers(
        {{ "Escuro (CartoDB Dark)": escuro, "Rico (CartoDB Voyager)": rico, "Satélite (Esri)": satelite, "Relevo/Topográfico": relevo, "OSM Padrão": osmPadrao }},
        overlaysMapa, {{ collapsed: false }}
    ).addTo(mapa);
    L.control.scale({{ metric: true, imperial: false }}).addTo(mapa);

    var destaqueMapa = L.circleMarker([0, 0], {{
        radius: 14, color: '{MARCA_ROXO}', weight: 3, fillOpacity: 0, opacity: 0,
    }}).addTo(mapa);

    // coluna estratigrafica ESTIMADA do ponto selecionado (furo virtual) --
    // blocos proporcionais a espessura, escala fixa (mesmo px/m pra todos os
    // pontos, da pra comparar a profundidade entre eles).
    var COLUNAS = {colunas_js};
    var COLUNA_ALTURA_PX = 460, COLUNA_MAX_M = {max_total_coluna:.0f};
    function atualizarColuna(id) {{
        var alvo = document.getElementById('coluna-estrat');
        var c = COLUNAS[id], r = DADOS_POR_ID[id];
        if (!c) {{
            alvo.innerHTML = '<div class="coluna-vazia">Sem coordenada para estimar a coluna deste item.</div>';
            return;
        }}
        var pxm = COLUNA_ALTURA_PX / COLUNA_MAX_M;
        var eixo = '', pilha = '';
        for (var d = 0; d <= COLUNA_MAX_M; d += 100) {{
            eixo += '<span style="top:' + (d * pxm) + 'px">' + d + '</span>';
        }}
        c.segs.forEach(function(s) {{
            var h = (s[2] - s[1]) * pxm;
            var rotulo = h >= 30 ? s[0] + '<br>' + s[1] + '–' + s[2] + ' m' : (h >= 15 ? s[0] : '');
            pilha += '<div class="coluna-seg" title="' + s[0] + ': ' + s[1] + '–' + s[2] + ' m (' + (s[2] - s[1]) + ' m)" '
                   + 'style="top:' + (s[1] * pxm) + 'px;height:' + h + 'px;background:' + s[3] + ';color:' + s[4] + '">'
                   + rotulo + '</div>';
        }});
        var lado = c.dique
            ? '<div class="coluna-dique" title="Dique (corpo subvertical, atravessa a coluna toda)" style="top:0;height:' + (c.total * pxm) + 'px"></div>'
            : '';
        var info = '<div class="coluna-titulo">' + (r ? r.nome : id) + '</div>'
                 + 'Cota do terreno: ' + c.z + ' m<br>Profundidade modelada: ' + c.total + ' m'
                 + (c.dique ? '<br><b>Dentro de polígono de dique</b> (barra verde ao lado)' : '')
                 + (c.d_sill !== null ? '<br>Sill mais próximo: ' + (c.d_sill / 1000).toFixed(2) + ' km' : '')
                 + (c.d_dique !== null ? '<br>Dique mais próximo: ' + (c.d_dique / 1000).toFixed(2) + ' km' : '')
                 + '<br><br><span style="opacity:0.7">Estimativa pelo modelo da seção 2D: plano de mergulho regional + erosão contra o relevo; espessuras da literatura (Teresina 350, Serra Alta 80, Irati 55, Palermo 100, Rio Bonito 269 m). Não é furo real.</span>';
        alvo.innerHTML = '<div class="coluna-corpo"><div class="coluna-eixo" style="height:' + COLUNA_ALTURA_PX + 'px">' + eixo + '</div>'
                       + '<div class="coluna-pilha" style="height:' + (c.total * pxm) + 'px">' + pilha + '</div>'
                       + '<div class="coluna-lado">' + lado + '</div>'
                       + '<div class="coluna-info">' + info + '</div></div>';
    }}

    function selecionarPorId(id, origemClique) {{
        atualizarColuna(id);
        document.querySelectorAll('.linha-dado.selecionada').forEach(function(el) {{ el.classList.remove('selecionada'); }});
        var linha = document.querySelector('.linha-dado[data-id="' + id + '"]');
        if (linha) {{
            linha.classList.add('selecionada');
            if (origemClique !== 'lista') linha.scrollIntoView({{block: 'nearest'}});
        }}
        atualizarDestaqueGraficos(id);
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
        if (r.geoq) {{
            linhas.push('<div class="linha-popup" style="margin-top:4px;"><b>Geoquímica (' + r.geoq.amostra + '):</b></div>');
            linhas.push('<div class="linha-popup">SiO₂ ' + r.geoq.SiO2.toFixed(1) + '% · TiO₂ ' + r.geoq.TiO2.toFixed(2)
                + '% · MgO ' + r.geoq.MgO.toFixed(2) + '% · P₂O₅ ' + r.geoq.P2O5.toFixed(2) + '%</div>');
            linhas.push('<div class="linha-popup">Sr ' + Math.round(r.geoq.Sr) + ' ppm · Zr ' + Math.round(r.geoq.Zr)
                + ' ppm · Y ' + r.geoq.Y.toFixed(1) + ' ppm</div>');
        }}
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
