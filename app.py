import numpy as np
import pandas as pd
from io import BytesIO
import matplotlib.pyplot as plt
import streamlit as st

st.set_page_config(page_title="Misión NASA 1962 – Método de Euler", page_icon="🚀", layout="wide")

# ---------------------------------------------------------------
# Paleta pastel en tonos azules (inspirada en el espacio)
# ---------------------------------------------------------------
AZUL_OSCURO = "#1B3A5C"
AZUL_ANALITICA = "#0F2A4A"
COLORES = {"Euler": "#5A9BE6", "Heun": "#2B5FAA", "RK4": "#17B1D6"}

plt.rcParams.update({
    "figure.facecolor": "#FFFFFF",
    "axes.facecolor": "#F5F9FF",
    "axes.edgecolor": "#BBD6F2",
    "axes.labelcolor": AZUL_OSCURO,
    "axes.titlecolor": AZUL_OSCURO,
    "xtick.color": AZUL_OSCURO,
    "ytick.color": AZUL_OSCURO,
    "grid.color": "#BBD6F2",
    "legend.facecolor": "#FFFFFF",
    "legend.edgecolor": "#BBD6F2",
})

def mostrar_tabla(df):
    """Muestra una DataFrame como tabla HTML con cabecera azul."""
    html = df.to_html(index=False, classes="tabla-azul", border=0,
                      float_format=lambda x: f"{x:.4f}").replace("\n", "")
    st.markdown(f'<div class="contenedor-tabla">{html}</div>', unsafe_allow_html=True)


st.markdown(
    """
    <style>
    .stApp {
        background: linear-gradient(160deg, #EAF4FF 0%, #F4F8FF 50%, #E3EEFB 100%);
        color: #1B3A5C;
    }
    [data-testid="stHeader"] { background: transparent; }
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #D6E8FA 0%, #E9F1FB 60%, #C9DFF5 100%);
    }
    [data-testid="stSidebar"] * { color: #1B3A5C; }
    h1 {
        color: #1B3A5C !important;
        border-bottom: 4px solid #6FA8E0;
        padding-bottom: 0.3rem;
    }
    h2, h3 { color: #24527F !important; }
    .stMarkdown, label, [data-testid="stCaptionContainer"] { color: #1B3A5C; }
    [data-testid="stMetric"] {
        background: rgba(255, 255, 255, 0.75);
        border: 1px solid #BBD6F2;
        border-left: 6px solid #4A90D9;
        border-radius: 12px;
        padding: 0.6rem 0.9rem;
    }
    .stTabs [data-baseweb="tab-list"] { gap: 6px; }
    .stTabs [data-baseweb="tab"] {
        background: #DCEBFA;
        border-radius: 10px 10px 0 0;
        padding: 0.4rem 1rem;
        color: #1B3A5C;
    }
    .stTabs [aria-selected="true"] {
        background: #C7DDF5 !important;
        color: #1B3A5C !important;
        font-weight: 600;
    }
    .stTabs [aria-selected="true"] p { color: #1B3A5C !important; }
    .stTabs [data-baseweb="tab-highlight"] { background-color: #8DB9E8 !important; }
    .stTabs [data-baseweb="tab-border"] { background-color: #D6E8FA !important; }
    [data-testid="stExpander"] {
        background: rgba(255, 255, 255, 0.65);
        border: 1px solid #BBD6F2;
        border-radius: 12px;
    }
    .stButton button, .stDownloadButton button {
        background: #4A90D9;
        color: #FFFFFF;
        border: 0;
        border-radius: 10px;
    }
    .stButton button:hover, .stDownloadButton button:hover {
        background: #3578BD;
        color: #FFFFFF;
    }
    /* Tablas con cabecera azul */
    .contenedor-tabla {
        max-height: 430px;
        overflow: auto;
        border: 1px solid #A9CBEE;
        border-radius: 10px;
        background: #FFFFFF;
        margin-bottom: 1rem;
    }
    .tabla-azul { border-collapse: collapse; width: 100%; font-size: 0.9rem; }
    .tabla-azul thead th {
        background: #4A90D9;
        color: #FFFFFF;
        padding: 8px 10px;
        text-align: right;
        position: sticky;
        top: 0;
    }
    .tabla-azul td {
        padding: 6px 10px;
        text-align: right;
        color: #1B3A5C;
        border-bottom: 1px solid #DCEBFA;
    }
    .tabla-azul tbody tr:nth-child(even) { background: #F0F7FF; }
    </style>
    """,
    unsafe_allow_html=True,
)

# ---------------------------------------------------------------
# Modelo:  m·dv/dt = T − m·g − k·v   →   dv/dt = a − b·v ,  dy/dt = v
# con a = T/m − g   y   b = k/m
# ---------------------------------------------------------------


def f(s, a, b):
    """Derivada del sistema. s = [v, y]."""
    return np.array([a - b * s[0], s[0]])


def paso_euler(s, h, a, b):
    return s + h * f(s, a, b)


def paso_heun(s, h, a, b):
    k1 = f(s, a, b)
    k2 = f(s + h * k1, a, b)
    return s + h / 2 * (k1 + k2)


def paso_rk4(s, h, a, b):
    k1 = f(s, a, b)
    k2 = f(s + h / 2 * k1, a, b)
    k3 = f(s + h / 2 * k2, a, b)
    k4 = f(s + h * k3, a, b)
    return s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


METODOS = {"Euler": paso_euler, "Heun": paso_heun, "RK4": paso_rk4}


def simular(metodo, h, tf, v0, y0, a, b):
    n = max(1, int(round(tf / h)))
    h = tf / n
    t = np.linspace(0, tf, n + 1)
    S = np.zeros((n + 1, 2))
    S[0] = [v0, y0]
    for i in range(n):
        S[i + 1] = METODOS[metodo](S[i], h, a, b)
    return t, S[:, 0], S[:, 1]


def exacta(t, v0, y0, a, b):
    t = np.asarray(t, dtype=float)
    if abs(b) < 1e-12:  # sin resistencia
        return v0 + a * t, y0 + v0 * t + a * t**2 / 2
    vinf = a / b
    e = np.exp(-b * t)
    v = vinf + (v0 - vinf) * e
    y = y0 + vinf * t + (v0 - vinf) * (1 - e) / b
    return v, y


# ---------------------------------------------------------------
# Interfaz
# ---------------------------------------------------------------
st.title("🚀 Misión espacial NASA 1962 — Despegue de un cohete")
st.caption("Métodos Numéricos · Euler, Heun y Runge-Kutta 4 frente a la solución analítica")


with st.expander("📘 Modelo matemático", expanded=False):
    st.markdown("Segunda ley de Newton con empuje, peso y resistencia lineal (sentido positivo hacia arriba):")

    st.latex(r"m\frac{dv}{dt}=T-mg-kv\;\Rightarrow\;\frac{dv}{dt}=\frac{T}{m}-g-\frac{k}{m}v,\qquad \frac{dy}{dt}=v")

    st.markdown("**Método de Euler:**")
    st.latex(r"v_{n+1}=v_n+h(a-bv_n),\qquad y_{n+1}=y_n+hv_n")

    st.markdown("**Solución analítica:**")
    st.latex(r"b \ne 0,\qquad a=\frac{T}{m}-g,\qquad b=\frac{k}{m}")
    st.latex(r"v(t)=\frac{a}{b}+(v_0-\frac{a}{b})e^{-bt}")
    st.latex(r"y(t)=y_0+\frac{a}{b}t+(v_0-\frac{a}{b})\frac{1-e^{-bt}}{b}")

    st.markdown(
        "Para los valores iniciales del enunciado "
        "(m = 1000 kg, T = 15000 N, g = 9.81 m/s², "
        "k = 100 kg/s, v₀ = 0 y y₀ = 0), se obtiene:"
    )
    st.latex(r"v(t)=51{,}9(1-e^{-0{,}1t})")
    st.latex(r"y(t)=51{,}9t-519(1-e^{-0{,}1t})")

    st.caption("Si la resistencia es nula (b = 0), se utiliza el caso especial de la solución analítica.")


# ----- Parámetros -----
st.sidebar.header("⚙️ Parámetros")
m = st.sidebar.number_input("Masa m (kg)", min_value=1.0, value=1000.0, step=50.0)
T = st.sidebar.number_input("Empuje T (N)", min_value=0.0, value=15000.0, step=500.0)
g = st.sidebar.number_input("Gravedad g (m/s²)", min_value=0.0, value=9.81, step=0.01)
k = st.sidebar.number_input("Resistencia lineal k (kg/s)", min_value=0.0, value=100.0, step=10.0)
v0 = st.sidebar.number_input("Velocidad inicial v(0) (m/s)", value=0.0)
y0 = st.sidebar.number_input("Altura inicial y(0) (m)", value=0.0)
tf = st.sidebar.number_input("Tiempo de simulación (s)", min_value=1.0, value=20.0, step=1.0)
h = st.sidebar.number_input("Paso h (s)", min_value=0.01, max_value=float(tf), value=2.0, step=0.25, format="%.2f")

st.sidebar.header("🎯 Requisitos de precisión")
tol_v = st.sidebar.number_input("Error máx. velocidad (m/s)", min_value=0.0, value=0.5)
tol_y = st.sidebar.number_input("Error máx. altura (m)", min_value=0.0, value=5.0)

st.sidebar.header("🧮 Métodos")
elegidos = st.sidebar.multiselect("Métodos a comparar", list(METODOS), default=list(METODOS))
if "Euler" not in elegidos:
    elegidos = ["Euler"] + elegidos

a = T / m - g
b = k / m
st.info(f"Ecuación resultante: **dv/dt = {a:.4f} − {b:.4f}·v**,  dy/dt = v   (h usado = {tf / max(1, round(tf / h)):.4f} s)")

# ----- Cálculos -----
res = {met: simular(met, h, tf, v0, y0, a, b) for met in elegidos}
v_ex_f, y_ex_f = exacta(tf, v0, y0, a, b)

tab1, tab2, tab3, tab4 = st.tabs(["📊 Resultados", "📋 Tabla de iteraciones", "📈 Gráficas", "🔬 Efecto del paso h"])

# ----- Resultados -----
with tab1:
    st.subheader(f"Resultados a los {tf:g} s")
    st.markdown(f"**Solución analítica:** v = {float(v_ex_f):.4f} m/s · y = {float(y_ex_f):.4f} m")
    cols = st.columns(len(elegidos))
    for col, met in zip(cols, elegidos):
        t, v, y = res[met]
        ev, ey = abs(float(v_ex_f) - v[-1]), abs(float(y_ex_f) - y[-1])
        with col:
            st.markdown(f"### {met}")
            st.metric("Velocidad (m/s)", f"{v[-1]:.4f}", f"error {ev:.4f}", delta_color="off")
            st.metric("Altura (m)", f"{y[-1]:.4f}", f"error {ey:.4f}", delta_color="off")
            if ev < tol_v and ey < tol_y:
                st.success("✅ Cumple la precisión")
            else:
                st.error("❌ No cumple: reduce el paso h")

    st.markdown("**Comparación de Euler con la solución analítica en instantes específicos**")
    t_e, v_e, y_e = res["Euler"]
    filas = []
    for tt in [2, 4, 10, 20]:
        idx = int(np.argmin(np.abs(t_e - tt)))
        if abs(t_e[idx] - tt) < 1e-9:
            ve, ye = exacta(tt, v0, y0, a, b)
            filas.append({"t (s)": tt, "v Euler": v_e[idx], "v exacta": float(ve), "y Euler": y_e[idx], "y exacta": float(ye)})
    if filas:
        mostrar_tabla(pd.DataFrame(filas).round(4))
    else:
        st.caption("Con este h no coinciden instantes de 2, 4, 10 y 20 s; prueba h = 2, 1, 0,5 o 0,25.")

# ----- Tabla -----
with tab2:
    met_t = st.selectbox("Método", elegidos)
    t, v, y = res[met_t]
    ve, ye = exacta(t, v0, y0, a, b)
    df = pd.DataFrame({
        "t (s)": t, f"v {met_t} (m/s)": v, "v exacta (m/s)": ve, "error v": np.abs(ve - v),
        f"y {met_t} (m)": y, "y exacta (m)": ye, "error y": np.abs(ye - y),
    }).round(4)
    mostrar_tabla(df)
        # Preparar archivo Excel
    excel_buffer = BytesIO()

    with pd.ExcelWriter(excel_buffer, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Resultados")

    # Botón para descargar Excel
    st.download_button(
        "⬇️ Descargar Excel",
        data=excel_buffer.getvalue(),
        file_name=f"tabla_{met_t}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )


# ----- Gráficas -----
with tab3:
    tt = np.linspace(0, tf, 400)
    vv, yy = exacta(tt, v0, y0, a, b)

    # Gráficas de velocidad y altura
    fig, ax = plt.subplots(1, 2, figsize=(7, 2.8), dpi=100)

    estilos = {
        "Euler": "o--",
        "Heun": "s--",
        "RK4": "^--"
    }

    for j, (nombre, ex, idx_) in enumerate([
        ("Velocidad (m/s)", vv, 1),
        ("Altura (m)", yy, 2)
    ]):
        ax[j].plot(tt, ex, "-", color=AZUL_ANALITICA, linewidth=2, label="Analítica")

        for met in elegidos:
            ax[j].plot(
                res[met][0],
                res[met][idx_],
                estilos[met],
                color=COLORES[met],
                markersize=3,
                linewidth=1,
                label=met
            )

        ax[j].set_xlabel("Tiempo (s)", fontsize=8)
        ax[j].set_ylabel(nombre, fontsize=8)
        ax[j].set_title(nombre + " vs tiempo", fontsize=9)
        ax[j].tick_params(axis="both", labelsize=7)
        ax[j].grid(alpha=0.3)
        ax[j].legend(fontsize=7)

    fig.tight_layout()
    st.pyplot(fig, use_container_width=False)
    plt.close(fig)

    # Gráficas de errores absolutos
    fig2, ax2 = plt.subplots(1, 2, figsize=(7, 2.4), dpi=100)

    for met in elegidos:
        t_, v_, y_ = res[met]
        ve_, ye_ = exacta(t_, v0, y0, a, b)

        ax2[0].plot(t_, np.abs(ve_ - v_), label=met, linewidth=1.5, color=COLORES[met])
        ax2[1].plot(t_, np.abs(ye_ - y_), label=met, linewidth=1.5, color=COLORES[met])

    titulos = [
        "Error absoluto en velocidad",
        "Error absoluto en altura"
    ]

    unidades = [
        "Error (m/s)",
        "Error (m)"
    ]

    for j, axx in enumerate(ax2):
        axx.set_title(titulos[j], fontsize=9)
        axx.set_xlabel("Tiempo (s)", fontsize=8)
        axx.set_ylabel(unidades[j], fontsize=8)
        axx.tick_params(axis="both", labelsize=7)
        axx.grid(alpha=0.3)
        axx.legend(fontsize=7)

    fig2.tight_layout()
    st.pyplot(fig2, use_container_width=False)
    plt.close(fig2)




# ----- Efecto del paso -----
with tab4:
    st.subheader("¿Qué pasa si reducimos h?")

    hs = sorted({2.0, 1.0, 0.5, 0.25, float(h)}, reverse=True)

    filas = []

    for met in elegidos:
        for hh in hs:
            # Calcular el número de intervalos y el paso efectivo
            n = max(1, int(round(tf / hh)))
            h_real = tf / n

            # Ejecutar la simulación
            _, v, y = simular(met, hh, tf, v0, y0, a, b)

            # Calcular errores absolutos
            error_v = abs(float(v_ex_f) - v[-1])
            error_y = abs(float(y_ex_f) - y[-1])

            filas.append({
                "Método": met,
                "h solicitado (s)": hh,
                "h efectivo (s)": h_real,
                "v final (m/s)": v[-1],
                "error v (m/s)": error_v,
                "y final (m)": y[-1],
                "error y (m)": error_y,
                "Cumple": "✅" if error_v < tol_v and error_y < tol_y else "❌"
            })

    # Mostrar tabla comparativa
    dfh = pd.DataFrame(filas)
    mostrar_tabla(dfh.round(4))

  
    # Gráfica de error en función del paso efectivo
    fig3, ax3 = plt.subplots(figsize=(4, 2.5), dpi=100)

    for met in elegidos:
        sub = dfh[dfh["Método"] == met].sort_values("h efectivo (s)")
        ax3.loglog(
            sub["h efectivo (s)"],
            sub["error y (m)"] + 1e-16,
            marker="o",
            color=COLORES[met],
            markersize=3,
            linewidth=1,
            label=met
        )

    ax3.set_xlabel("Paso efectivo h (s)", fontsize=7)
    ax3.set_ylabel("Error en altura (m)", fontsize=7)
    ax3.set_title("Error en altura vs paso de integración", fontsize=8)

    # Tamaño uniforme de los números de ambos ejes
    ax3.tick_params(
        axis="both",
        which="both",
        labelsize=5
    )

    for etiqueta in ax3.xaxis.get_majorticklabels():
        etiqueta.set_fontsize(5)

    for etiqueta in ax3.yaxis.get_majorticklabels():
        etiqueta.set_fontsize(5)

    ax3.grid(alpha=0.25, which="both", linewidth=0.5)
    ax3.legend(fontsize=6, markerscale=0.7)

    fig3.tight_layout()
    st.pyplot(fig3, use_container_width=False)


    st.markdown(
        """
        **¿Cómo interpretar los resultados?**

        - **Euler:** es un método de orden 1. Al reducir el paso a la mitad,
          el error global tiende aproximadamente a reducirse a la mitad.
        - **Heun:** es de orden 2. El error global tiende aproximadamente
          a reducirse a una cuarta parte.
        - **RK4:** es de orden 4. El error global tiende aproximadamente
          a reducirse a una dieciseisava parte.

        Para un mismo paso de integración, RK4 generalmente proporciona
        una aproximación más precisa que Heun y Euler.
        """
    )
