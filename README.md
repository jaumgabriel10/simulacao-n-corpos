# Simulação Gravitacional de N Corpos

Projeto de Física Computacional dedicado à simulação numérica da dinâmica gravitacional de um sistema de N corpos, utilizando métodos de integração temporal e aproximações numéricas para o cálculo das interações gravitacionais.

## Sobre o projeto

O objetivo é estudar o comportamento de sistemas gravitacionais compostos por múltiplos corpos, observando fenômenos como formação de estruturas, interações gravitacionais, estabilidade orbital e ejeção de partículas.

As condições iniciais são geradas aleatoriamente, permitindo investigar diferentes configurações e comportamentos dinâmicos do sistema.

## Métodos numéricos

### 1. Gravitação newtoniana

A interação gravitacional entre os corpos é descrita pela Lei da Gravitação Universal de Newton:

\[
\vec{F}_{ij} = G\frac{m_i m_j}{r_{ij}^3}\vec{r}_{ij}
\]

É utilizado um parâmetro de suavização (*softening*) para evitar divergências numéricas em distâncias muito pequenas.

### 2. Integração temporal — Velocity Verlet

A evolução temporal do sistema é calculada pelo método Velocity Verlet, um integrador simplético frequentemente utilizado em simulações de sistemas físicos conservativos.

### 3. Algoritmo Barnes–Hut

Para otimizar o cálculo das forças gravitacionais, é implementado o algoritmo Barnes–Hut, baseado em uma estrutura de árvore quadtree.

Em vez de calcular individualmente a interação entre todos os corpos, o algoritmo aproxima grupos suficientemente distantes por uma única massa equivalente localizada em seu centro de massa.

O critério de abertura utilizado é:

\[
\frac{s}{d}<\theta
\]

Em que:

- \(s\): tamanho da região representada pelo nó;
- \(d\): distância entre o corpo e o centro de massa do nó;
- \(\theta\): parâmetro de abertura.

Essa aproximação reduz o custo computacional em comparação com o cálculo direto de todas as interações.

## Funcionalidades

- Geração aleatória de sistemas com N corpos;
- Interação gravitacional entre partículas;
- Integração temporal pelo método Velocity Verlet;
- Cálculo de forças utilizando Barnes–Hut;
- Visualização dinâmica das partículas;
- Acompanhamento do centro de massa do sistema;
- Identificação e contagem de partículas ejetadas;
- Visualização do potencial gravitacional em uma malha;
- Monitoramento da energia cinética, potencial e total;
- Monitoramento do momento angular;
- Análise dos erros relativos de conservação.

## Parâmetros da simulação

| Parâmetro | Valor |
|---|---:|
| Constante gravitacional \(G\) | 1 |
| Suavização gravitacional \(\epsilon\) | 0,08 |
| Passo temporal \(\Delta t\) | 0,005 |
| Parâmetro Barnes–Hut \(\theta\) | 0,6 |
| Profundidade máxima da árvore | 20 |
| Raio inicial de distribuição | 5 |

As grandezas são expressas em unidades normalizadas.

## Tecnologias utilizadas

- **Python** — implementação da simulação;
- **NumPy** — operações numéricas e vetoriais;
- **Matplotlib** — visualização e interface gráfica.

## Relatório

O relatório em PDF apresenta a fundamentação física, os métodos numéricos empregados, os parâmetros adotados, as funcionalidades implementadas e uma discussão das limitações da simulação.

Consulte o arquivo `Relatorio_ncorpos_joaogabriel.pdf` neste repositório.

## Objetivo acadêmico

Este projeto foi desenvolvido no contexto de estudos em Física Computacional, com o propósito de explorar métodos numéricos aplicados a sistemas gravitacionais de muitos corpos e compreender os desafios computacionais envolvidos na evolução dinâmica desses sistemas.

---

**Autor:** João Gabriel Fernandes Romero

**Área:** Física Computacional

**Instituição:** Universidade Federal do Rio Grande — FURG
