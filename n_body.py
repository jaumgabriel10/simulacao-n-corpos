import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, Slider

#como entrar no ambiente:
#  cd "E:\NeymarJr\Faculdade\Física Computacional"
# Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
#.\nbody312\Scripts\Activate.ps1
#como rodar o arquivo: python n_body.py

#Obs.: NÃO MEXA NA ANIMAÇÃO, CAMERA, EJEÇÃO e ÁRVORE. A IA ME AJUDOU COM ISSO, E EU NÃO QUERO QUEBRAR NADA. SÓ MEXA NOS PARÂMETROS ABAIXO.
#obs.: Ate 50 corpos da pra ter uma ideia, passou disso, a simulação fica lenta. Use até no max 100

G = 1.0
dt = 0.005
passos_por_frame = 20
intervalo = 15
suavizacao = 0.08
theta = 0.6
raio_ejecao = 12.0
limite_tela = 8.0
atualizar_campo_a_cada = 5
rng = np.random.default_rng()

while True:
    try:
        N_inicial = int(input("Quantos corpos deseja simular? "))
        if N_inicial > 0:
            break
        print("Digite um número maior que zero.")
    except ValueError:
        print("Digite um número inteiro válido.")

#Obs.: essa parte foi inteiramente com ajuda de IA - Aproximação de Barnes-Hut
class No:
    def __init__(self, indices, centro, meia_largura, profundidade=0):
        self.indices = indices
        self.centro = centro
        self.meia_largura = meia_largura
        self.profundidade = profundidade
        self.filhos = []
        self.massa = np.sum(massas[indices])
        self.com = np.sum(posicoes[indices] * massas[indices, None], axis=0) / self.massa

        if len(indices) > 1 and profundidade < 20:
            self.subdividir()

    def subdividir(self):
        x, y = posicoes[self.indices].T
        cx, cy = self.centro
        quadrantes = [
            self.indices[(x < cx) & (y < cy)],
            self.indices[(x >= cx) & (y < cy)],
            self.indices[(x < cx) & (y >= cy)],
            self.indices[(x >= cx) & (y >= cy)]
        ]

        if np.max(np.ptp(posicoes[self.indices], axis=0)) < 1e-12:
            return

        deslocamento = self.meia_largura / 2
        centros = [
            (cx - deslocamento, cy - deslocamento),
            (cx + deslocamento, cy - deslocamento),
            (cx - deslocamento, cy + deslocamento),
            (cx + deslocamento, cy + deslocamento)
        ]

        for indices, centro in zip(quadrantes, centros):
            if len(indices):
                self.filhos.append(
                    No(indices, centro, deslocamento, self.profundidade + 1)
                )

#Obs.: essa parte foi inteiramente com ajuda de IA - Aproximação de Barnes-Hut
def construir_arvore():
    minimo = np.min(posicoes, axis=0)
    maximo = np.max(posicoes, axis=0)
    centro = (minimo + maximo) / 2
    meia_largura = max(np.max(maximo - minimo) / 2, 1.0) * 1.01
    return No(np.arange(len(massas)), centro, meia_largura)


def contem_corpo(no, indice):
    deslocamento = np.abs(posicoes[indice] - no.centro)
    return np.all(deslocamento <= no.meia_largura)


def forca_no_corpo(no, indice):
    diferenca = no.com - posicoes[indice]
    distancia2 = np.dot(diferenca, diferenca) + suavizacao**2

    if not no.filhos:
        resultado = np.zeros(2)
        for j in no.indices:
            if j == indice:
                continue
            d = posicoes[j] - posicoes[indice]
            r2 = np.dot(d, d) + suavizacao**2
            resultado += G * massas[j] * d / r2**1.5
        return resultado

    distancia = np.sqrt(distancia2)
    tamanho = 2 * no.meia_largura

    if not contem_corpo(no, indice) and tamanho / distancia < theta:
        return G * no.massa * diferenca / distancia2**1.5

    resultado = np.zeros(2)
    for filho in no.filhos:
        resultado += forca_no_corpo(filho, indice)
    return resultado


def calcular_aceleracoes():
    if len(massas) == 0:
        return np.empty((0, 2))

    arvore = construir_arvore()
    aceleracoes_calc = np.zeros_like(posicoes)

    for i in range(len(massas)):
        aceleracoes_calc[i] = forca_no_corpo(arvore, i)

    return aceleracoes_calc


def potencial_no_corpo(no, indice):
    diferenca = no.com - posicoes[indice]
    distancia2 = np.dot(diferenca, diferenca) + suavizacao**2

    if not no.filhos:
        potencial = 0.0
        for j in no.indices:
            if j == indice:
                continue
            d = posicoes[j] - posicoes[indice]
            r = np.sqrt(np.dot(d, d) + suavizacao**2)
            potencial -= G * massas[j] / r
        return potencial

    distancia = np.sqrt(distancia2)
    tamanho = 2 * no.meia_largura

    if not contem_corpo(no, indice) and tamanho / distancia < theta:
        return -G * no.massa / distancia

    return sum(potencial_no_corpo(filho, indice) for filho in no.filhos)


def calcular_energias(pos, vel):
    K = 0.5 * np.sum(massas[:, None] * vel**2)

    if len(massas) < 2:
        return K, 0.0, K

    arvore = construir_arvore()
    potenciais = np.array([
        potencial_no_corpo(arvore, i) for i in range(len(massas))
    ])
    U = 0.5 * np.sum(massas * potenciais)

    return K, U, K + U


def calcular_momento_angular(pos, vel):
    if len(massas) == 0:
        return 0.0
    return np.sum(
        massas * (pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0])
    )


def gerar_condicoes_iniciais(n):
    m = rng.uniform(0.5, 2.0, n)

    angulos = rng.uniform(0, 2 * np.pi, n)
    raios = 5 * np.sqrt(rng.random(n))
    p = np.column_stack((raios * np.cos(angulos), raios * np.sin(angulos)))

    direcoes = rng.uniform(0, 2 * np.pi, n)
    intensidades = rng.uniform(0.1, 1.2, n)
    v = np.column_stack((
        intensidades * np.cos(direcoes),
        intensidades * np.sin(direcoes)
    ))

    identificadores = [set([i]) for i in range(n)]
    return m, p, v, identificadores


def detectar_ejecoes():
    global ejetados

    if len(massas) < 2:
        return

    massa_total = np.sum(massas)
    centro_massa = np.sum(posicoes * massas[:, None], axis=0) / massa_total
    velocidade_cm = np.sum(velocidades * massas[:, None], axis=0) / massa_total

    arvore = construir_arvore()

    for i in range(len(massas)):
        r = posicoes[i] - centro_massa
        v = velocidades[i] - velocidade_cm
        distancia = np.linalg.norm(r)

        if distancia <= raio_ejecao or np.dot(r, v) <= 0:
            continue

        energia_especifica = 0.5 * np.dot(v, v) + potencial_no_corpo(arvore, i)

        if energia_especifica > 0:
            ejetados.update(identificadores[i])


def calcular_centro_massa():
    massa_total = np.sum(massas)
    if massa_total == 0:
        return np.zeros(2)
    return np.sum(posicoes * massas[:, None], axis=0) / massa_total


def calcular_campo(centro=None):
    if centro is None:
        centro = calcular_centro_massa()

    x = np.linspace(centro[0] - limite_tela, centro[0] + limite_tela, 65)
    y = np.linspace(centro[1] - limite_tela, centro[1] + limite_tela, 65)
    X, Y = np.meshgrid(x, y)
    potencial = np.zeros_like(X)

    for m, p in zip(massas, posicoes):
        r2 = (X - p[0])**2 + (Y - p[1])**2 + suavizacao**2
        potencial -= G * m / np.sqrt(r2)

    return potencial


def atualizar_referencias():
    global K0, U0, E0, L0, escala_E, escala_L
    K0, U0, E0 = calcular_energias(posicoes, velocidades)
    L0 = calcular_momento_angular(posicoes, velocidades)
    escala_E = max(abs(E0), 1e-12)
    escala_L = max(abs(L0), 1e-12)


massas, posicoes, velocidades, identificadores = gerar_condicoes_iniciais(N_inicial)
ejetados = set()
aceleracoes = calcular_aceleracoes()
atualizar_referencias()

tempo = 0.0
pausado = False
campo_ativo = True
historico_tempo = [0.0]
historico_erro_tempo = [0.0]
historico_K = []
historico_U = []
historico_E = []
historico_erro_E = []
historico_erro_L = []

#eu tinha esquecido de registrar, rs
def registrar_historico():
    K, U, E = calcular_energias(posicoes, velocidades)
    L = calcular_momento_angular(posicoes, velocidades)

    historico_tempo.append(tempo)
    historico_K.append(K)
    historico_U.append(U)
    historico_E.append(E)
    historico_erro_E.append(max(abs((E - E0) / escala_E), 1e-15))
    historico_erro_L.append(max(abs((L - L0) / escala_L), 1e-15))


fig = plt.figure(figsize=(14, 8))
ax_orbita = fig.add_subplot(2, 2, (1, 3))
ax_energia = fig.add_subplot(2, 2, 2)
ax_erros = fig.add_subplot(2, 2, 4)

fig.subplots_adjust(bottom=0.19, top=0.91, wspace=0.28, hspace=0.32)

ax_orbita.set_title("Movimento dos corpos")
ax_orbita.set_xlabel("Posição X")
ax_orbita.set_ylabel("Posição Y")
ax_orbita.set_xlim(-limite_tela, limite_tela)
ax_orbita.set_ylim(-limite_tela, limite_tela)
ax_orbita.set_aspect("equal")
ax_orbita.grid(alpha=0.25)

centro_inicial = calcular_centro_massa()
ax_orbita.set_xlim(-limite_tela, limite_tela)
ax_orbita.set_ylim(-limite_tela, limite_tela)

campo = ax_orbita.imshow(
    calcular_campo(centro_inicial),
    extent=(-limite_tela, limite_tela, -limite_tela, limite_tela),
    origin="lower",
    cmap="magma",
    alpha=0.42,
    interpolation="bilinear",
    zorder=0
)

velocidades_modulo = np.linalg.norm(velocidades, axis=1)
pontos = ax_orbita.scatter(
    posicoes[:, 0] - centro_inicial[0],
    posicoes[:, 1] - centro_inicial[1],
    s=25 + 15 * np.sqrt(massas),
    c=velocidades_modulo,
    cmap="turbo",
    vmin=0,
    vmax=max(np.max(velocidades_modulo), 1.0),
    edgecolors="white",
    linewidths=0.35,
    zorder=2
)

barra_cores = fig.colorbar(pontos, ax=ax_orbita, fraction=0.046, pad=0.04)
barra_cores.set_label("Velocidade")

linha_K, = ax_energia.plot([], [], label="Cinética")
linha_U, = ax_energia.plot([], [], label="Potencial")
linha_E, = ax_energia.plot([], [], label="Total")
ax_energia.set_title("Energia do sistema")
ax_energia.set_xlabel("Tempo")
ax_energia.set_ylabel("Energia")
ax_energia.legend()
ax_energia.grid(alpha=0.3)

linha_erro_E, = ax_erros.plot([], [], label="Energia")
linha_erro_L, = ax_erros.plot([], [], label="Momento angular")
ax_erros.set_title("Erro relativo em relação ao início")
ax_erros.set_xlabel("Tempo")
ax_erros.set_ylabel("Erro relativo")
ax_erros.set_yscale("log")
ax_erros.set_ylim(1e-15, 1e2)
ax_erros.legend()
ax_erros.grid(alpha=0.3)

ax_pausa = fig.add_axes([0.08, 0.055, 0.12, 0.055])
ax_reiniciar = fig.add_axes([0.22, 0.055, 0.12, 0.055])
ax_campo = fig.add_axes([0.36, 0.055, 0.13, 0.055])
ax_velocidade = fig.add_axes([0.58, 0.065, 0.30, 0.035])

botao_pausa = Button(ax_pausa, "Pausar")
botao_reiniciar = Button(ax_reiniciar, "Reiniciar")
botao_campo = Button(ax_campo, "Campo: ON")
slider_velocidade = Slider(
    ax_velocidade, "Passos/frame", 1, 100,
    valinit=passos_por_frame, valstep=1
)


def integrar():
    global posicoes, velocidades, aceleracoes, tempo

    novas_posicoes = (
        posicoes + velocidades * dt + 0.5 * aceleracoes * dt**2
    )

    posicoes = novas_posicoes
    novas_aceleracoes = calcular_aceleracoes()

    novas_velocidades = (
        velocidades + 0.5 * (aceleracoes + novas_aceleracoes) * dt
    )

    velocidades = novas_velocidades
    aceleracoes = novas_aceleracoes
    tempo += dt


#IA na animação, rs
def alternar_pausa(evento):
    global pausado

    if pausado:
        animacao.event_source.start()
        botao_pausa.label.set_text("Pausar")
    else:
        animacao.event_source.stop()
        botao_pausa.label.set_text("Continuar")

    pausado = not pausado
    fig.canvas.draw_idle()


def alterar_velocidade(valor):
    global passos_por_frame
    passos_por_frame = int(valor)


def alternar_campo(evento):
    global campo_ativo
    campo_ativo = not campo_ativo
    campo.set_visible(campo_ativo)
    botao_campo.label.set_text("Campo: ON" if campo_ativo else "Campo: OFF")
    fig.canvas.draw_idle()


def reiniciar(evento):
    global massas, posicoes, velocidades, aceleracoes, identificadores
    global tempo, ejetados
    global historico_tempo, historico_erro_tempo
    global historico_K, historico_U, historico_E
    global historico_erro_E, historico_erro_L

    massas, posicoes, velocidades, identificadores = gerar_condicoes_iniciais(N_inicial)
    aceleracoes = calcular_aceleracoes()
    tempo = 0.0
    ejetados = set()
    
    atualizar_referencias()

    historico_tempo = [0.0]
    historico_erro_tempo = [0.0]
    historico_K = []
    historico_U = []
    historico_E = []
    historico_erro_E = []
    historico_erro_L = []

    centro = atualizar_camera()
    atualizar_pontos()
    campo.set_data(calcular_campo(centro))
    campo.set_extent((-limite_tela, limite_tela, -limite_tela, limite_tela))
    linha_K.set_data([], [])
    linha_U.set_data([], [])
    linha_E.set_data([], [])
    linha_erro_E.set_data([], [])
    linha_erro_L.set_data([], [])
    ax_energia.set_xlim(0, 1)
    ax_erros.set_xlim(0, 1)
    fig.canvas.draw_idle()

#IA me ajudou com o negocio da camera
def atualizar_camera():
    centro = calcular_centro_massa()
    largura = ax_orbita.get_xlim()[1] - ax_orbita.get_xlim()[0]
    altura = ax_orbita.get_ylim()[1] - ax_orbita.get_ylim()[0]

    #CM no centro da tela
    ax_orbita.set_xlim(-largura / 2, largura / 2)
    ax_orbita.set_ylim(-altura / 2, altura / 2)
    return centro


def atualizar_pontos():
    centro = calcular_centro_massa()
    pontos.set_offsets((posicoes - centro).reshape(-1, 2))
    pontos.set_sizes(25 + 15 * np.sqrt(massas))
    velocidades_modulo = np.linalg.norm(velocidades, axis=1)
    pontos.set_array(velocidades_modulo)

    vmax = max(np.max(velocidades_modulo) if len(velocidades_modulo) else 1, 1.0)
    pontos.set_clim(0, vmax)


def zoom(evento):
    if evento.inaxes != ax_orbita or evento.xdata is None or evento.ydata is None:
        return

    fator = 0.85 if evento.button == "up" else 1.18
    xlim = ax_orbita.get_xlim()
    ylim = ax_orbita.get_ylim()

    largura = (xlim[1] - xlim[0]) * fator
    altura = (ylim[1] - ylim[0]) * fator

    fx = (evento.xdata - xlim[0]) / (xlim[1] - xlim[0])
    fy = (evento.ydata - ylim[0]) / (ylim[1] - ylim[0])

    ax_orbita.set_xlim(
        evento.xdata - largura * fx,
        evento.xdata + largura * (1 - fx)
    )
    ax_orbita.set_ylim(
        evento.ydata - altura * fy,
        evento.ydata + altura * (1 - fy)
    )
    fig.canvas.draw_idle()


def atualizar(frame):
    global tempo

    for _ in range(passos_por_frame):
        integrar()

    if frame % 10 == 0:
        detectar_ejecoes()

    centro = atualizar_camera()
    atualizar_pontos()

    if campo_ativo and frame % atualizar_campo_a_cada == 0:
        campo.set_data(calcular_campo(centro))
        campo.set_extent((-limite_tela, limite_tela, -limite_tela, limite_tela))

    if frame % 5 == 0:
        K, U, E = calcular_energias(posicoes, velocidades)
        L = calcular_momento_angular(posicoes, velocidades)

        historico_tempo.append(tempo)
        historico_K.append(K)
        historico_U.append(U)
        historico_E.append(E)
        historico_erro_tempo.append(tempo)
        historico_erro_E.append(max(abs((E - E0) / escala_E), 1e-15))
        historico_erro_L.append(max(abs((L - L0) / escala_L), 1e-15))

        linha_K.set_data(historico_tempo[1:], historico_K)
        linha_U.set_data(historico_tempo[1:], historico_U)
        linha_E.set_data(historico_tempo[1:], historico_E)
        linha_erro_E.set_data(historico_erro_tempo[1:], historico_erro_E)
        linha_erro_L.set_data(historico_erro_tempo[1:], historico_erro_L)

        ax_energia.set_xlim(0, max(tempo, dt))
        ax_erros.set_xlim(0, max(tempo, dt))
        ax_energia.relim()
        ax_energia.autoscale_view(scalex=False, scaley=True)
        ax_erros.relim()
        ax_erros.autoscale_view(scalex=False, scaley=True)

    fig.suptitle(
        f"Simulação gravitacional | Corpos: {len(massas)} | "
        f"Ejetados: {len(ejetados)} | "
        f"Tempo: {tempo:.2f}"
    )

    return (
        pontos, campo, linha_K, linha_U, linha_E,
        linha_erro_E, linha_erro_L
    )


botao_pausa.on_clicked(alternar_pausa)
botao_reiniciar.on_clicked(reiniciar)
botao_campo.on_clicked(alternar_campo)
slider_velocidade.on_changed(alterar_velocidade)
fig.canvas.mpl_connect("scroll_event", zoom)

animacao = FuncAnimation(
    fig,
    atualizar,
    interval=intervalo,
    blit=False,
    cache_frame_data=False
)

plt.show()

Kf, Uf, Ef = calcular_energias(posicoes, velocidades)
Lf = calcular_momento_angular(posicoes, velocidades)

print("\nResultados finais")
print(f"Tempo simulado: {tempo:.5f}")
print(f"Corpos restantes: {len(massas)}")
print(f"Corpos ejetados: {len(ejetados)}")
print(f"Energia final: {Ef:.10e}")
print(f"Momento angular final: {Lf:.10e}")
