import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

st.set_page_config(page_title="Misión NASA 1962 – Método de Euler", page_icon="🚀", layout="wide")

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
    st.latex(r"m\frac{dv}{dt}=T-mg-kv \;\Rightarrow\; \frac{dv}{dt}=\frac{T}{m}-g-\frac{k}{m}v,\qquad \frac{dy}{dt}=v")
    st.markdown("**Euler:**")
    st.latex(r"v_{n+1}=v_n+h\,(a-bv_n),\qquad y_{n+1}=y_n+h\,v_n")
    st.markdown("**Solución analítica** (con $a=T/m-g$, $b=k/m$):")
    st.latex(r"v(t)=\frac{a}{b}\left(1-e^{-bt}\right)+v_0e^{-bt},\qquad y(t)=y_0+\frac{a}{b}t+\left(v_0-\frac{a}{b}\right)\frac{1-e^{-bt}}{b}")
    st.markdown("Con los datos del enunciado: $v(t)=51{,}9\\,(1-e^{-0{,}1t})$ y $y(t)=51{,}9\\,t-519\\,(1-e^{-0{,}1t})$.")

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

    st.markdown("**Velocidad y altura en instantes específicos (Euler)**")
    t_e, v_e, y_e = res["Euler"]
    filas = []
    for tt in [2, 4, 10, 20]:
        idx = int(np.argmin(np.abs(t_e - tt)))
        if abs(t_e[idx] - tt) < 1e-9:
            ve, ye = exacta(tt, v0, y0, a, b)
            filas.append({"t (s)": tt, "v Euler": v_e[idx], "v exacta": float(ve), "y Euler": y_e[idx], "y exacta": float(ye)})
    if filas:
        st.dataframe(pd.DataFrame(filas).round(4), hide_index=True)
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
    st.dataframe(df, hide_index=True)
    st.download_button("⬇️ Descargar CSV", df.to_csv(index=False).encode("utf-8"), f"tabla_{met_t}.csv", "text/csv")

# ----- Gráficas -----
with tab3:
    tt = np.linspace(0, tf, 400)
    vv, yy = exacta(tt, v0, y0, a, b)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    estilos = {"Euler": "o--", "Heun": "s--", "RK4": "^--"}
    for j, (nombre, ex, idx_) in enumerate([("Velocidad (m/s)", vv, 1), ("Altura (m)", yy, 2)]):
        ax[j].plot(tt, ex, "k-", lw=2, label="Analítica")
        for met in elegidos:
            ax[j].plot(res[met][0], res[met][idx_], estilos[met], ms=4, label=met)
        ax[j].set_xlabel("t (s)")
        ax[j].set_ylabel(nombre)
        ax[j].set_title(nombre + " vs tiempo")
        ax[j].grid(alpha=0.3)
        ax[j].legend()
    st.pyplot(fig)

    fig2, ax2 = plt.subplots(1, 2, figsize=(12, 3.8))
    for met in elegidos:
        t_, v_, y_ = res[met]
        ve_, ye_ = exacta(t_, v0, y0, a, b)
        ax2[0].plot(t_, np.abs(ve_ - v_), label=met)
        ax2[1].plot(t_, np.abs(ye_ - y_), label=met)
    for axx, tl in zip(ax2, ["Error absoluto en velocidad (m/s)", "Error absoluto en altura (m)"]):
        axx.set_title(tl)
        axx.set_xlabel("t (s)")
        axx.grid(alpha=0.3)
        axx.legend()
    st.pyplot(fig2)


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
    st.dataframe(dfh.round(4), hide_index=True, use_container_width=True)

    # Gráfica de error en función del paso efectivo
    fig3, ax3 = plt.subplots(figsize=(2.5,1.5), dpi=100)
    
    for met in elegidos:
        sub = dfh[dfh["Método"] == met].sort_values("h efectivo (s)")

        ax3.loglog(
            sub["h efectivo (s)"],
            sub["error y (m)"] + 1e-16,
            "o-",
            label=met
        )

    ax3.set_xlabel("Paso efectivo h (s)")
    ax3.set_ylabel("Error en altura (m)")
    ax3.set_title("Error en altura vs paso de integración")
    ax3.grid(alpha=0.3, which="both")
    ax3.legend()

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

