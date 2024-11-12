# Bibliotecas para trabalhos numéricos
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

# Parâmetros conforme a numeração de grupo: 15B
T_inf   = 25.0      # Temperatura ao longe, em ºC
h       = 100.0     # Coef. de transf. de calor convectivo, em W/m²K
q0      = 8e4       # Fluxo de calor injetado, em W/m²

# Parâmetros gerais:
t_h     = 10        # Tempo de aquecimento, em s
w1      = 44e-3     # Comprimento da fita, em m
w2      = 132e-3    # Comprimento da placa, em m
d       = 1.25e-3   # Espessura da placa, em m
k       = 60.0      # Condutividade da fita, em W/mK
rho     = 7850.0    # Densidade da fita, em kg/m³
Cp      = 435.0     # Calor específico da fita, em J/kgK

# Parâmetros auxiliares
L       = w2/2                  # Metade do comprimento da placa (problema simétrico)
Bi      = (h*L**2)/(k*d)        # Número de Biot
alpha   = k/(rho*Cp)            # Difusividade térmica da fita
p       = w2/w1                 # Comprimento relativo entre a fita e a placa
q0_h    = (q0*L**2)/(k*d*T_inf) # Expressão do calor introduzido pela fita

# Parâmetros de integração
RES     = 100   # Nº de elementos para integração (resolução)
t_f     = 40    # Tempo de integração, em s
t       = np.linspace(0,t_f,RES) # Passo de tempo, em s
x       = np.linspace(0,L,RES) # Passo de distância, em m

# Parâmetros de solução adimensionais
t_      = (4*alpha*t)/(w2**2)   # Tempo adimensional
x_      = (2*x)/w2              # Distância adimensional

# Definição da matriz A conforme já resolvido
A = np.zeros([RES,RES])
A[0][0] = -2/3
A[0][1] =  2/3
for i in range(1,RES-1):
    A[i][i-1] =  1
    A[i][i]   = -2
    A[i][i+1] =  1
A[-1][-2] =  1
A[-1][-1] = -2

A = A/(2/RES)**2

# Função composta do calor introduzido
# Matriz q_
q_ = np.array([[(q0_h if abs(xi) <= 1/p and ti <= (4*alpha*t_h)/(w2**2) else 0) for xi in x_] for ti in t_])

# Implicitamente...
I = np.identity(RES)
inv = np.linalg.inv((1+2*Bi*(1/RES))*I-(1/RES)*A)

theta = np.zeros([RES,RES])

for i in range(RES-1):
    theta[i+1] = np.matmul(inv,(q_[i]*(1/RES) + theta[i]))


T = T_inf*(theta+1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x,t)

ax.set_xlabel('Distância x (mm)')
ax.set_ylabel('Tempo (s)')
ax.set_zlabel('Temperatura (ºC)')

my_cmap = plt.get_cmap('turbo')

ax.plot_surface(X*1000, Y, T, cmap=my_cmap,
                       edgecolor ='none')

ax.view_init(elev=30, azim=30, roll=0)
plt.show()
