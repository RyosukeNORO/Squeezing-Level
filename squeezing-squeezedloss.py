####スクイーズド光の損失がどの程度スクイージングレベルに影響を与えるかというプログラム（励起光損失は0と仮定している）
###スクイーズド光の損失があるときのスクイーズンぐレベルを計算するプログラム (励起光の損失は0としている)###
import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt

#定数
c = 299792458E6                 # 光速(um/s)
eps = 8.85418782E-18            # 真空の誘電率(F/um)

#変数
div = 1001
z = np.linspace(0, 40000, div)  # 素子長(um)
I = 0.1                        # ポンプ光パワー(W)
A = 10                          # 導波路断面積(um^2)
d33 = (2/np.pi)*13.8            # 非線形光学定数(pm/V)
d = d33*10E-6                   # 非線形光学定数(um/V)
npu = 2.271213253               # ポンプ光の感じる屈折率
ns = 2.147041829                # シグナル光の感じる屈折率
lamp = 0.405                    # ポンプ光波長(um)
lams = lamp*2                   # シグナル光波長(um)
kp = 2*np.pi*npu/lamp           # ポンプ光波数(1/um)
ks = 2*np.pi*ns/lams            # シグナル光波数(1/um)
omes = ks*c/ns                  # シグナル光の角周波数(rad/s)
kappa = 3.63                     # 非線形接合係数kappa (W^(1/2)cm^(-1))
kap = kappa / 10000
a = [0.23,0.35,0.7,1.5,3.5]           # 損失(dB/cm)
print(str(I*100/A)+'MW/cm^2')

# gの計算
#g = (2*np.pi*2*d)/npu**(3/2) * np.sqrt((8*np.pi*I)/(c*A)) * (ks/npu)
g = kap * np.sqrt(I) # np.sqrt((2*omes**2*d**2*I)/(ns**2*npu*eps*c**3*A)) #パラメトリックゲイン
print("g = {} /um".format(g))

#-----------------------------------------------------------------------------
# 0
gaa = a[0] # 導波路損失(dB/cm)
gama = np.log(10**(gaa*1E-5)) # 導波路損失(/um)
Sa = [0 for a in range(div)]
Sfa = [0 for a in range(div)]
Smaxa = 0

i = 0
while i < div:
    Sa[i] = 1 - np.exp(-gama*z[i]) + np.exp(-gama*z[i])*np.exp(-2*g*z[i])#(gama + 2*g*np.exp(-(gama+2*g)*z[i]))/(gama+2*g)
    Sfa[i] = -10*np.log10(Sa[i])
    if Smaxa < Sfa[i]:
        Smaxa = Sfa[i]
        zmaxa = z[i]
    i += 1

#-----------------------------------------------------------------------------
# 1
gab = a[1] # 導波路損失(dB/cm)
gamb = np.log(10**(gab*1E-5)) # 導波路損失(/um)
Sb = [0 for a in range(div)]
Sfb = [0 for a in range(div)]
Smaxb = 0

i = 0
while i < div:
    Sb[i] = 1 - np.exp(-gamb*z[i]) + np.exp(-gamb*z[i])*np.exp(-2*g*z[i])#(gamb + 2*g*np.exp(-(gamb+2*g)*z[i]))/(gamb+2*g)
    Sfb[i] = -10*np.log10(Sb[i])
    if Smaxb < Sfb[i]:
        Smaxb = Sfb[i]
        zmaxb = z[i]
    i += 1

#-----------------------------------------------------------------------------
# 2
gac = a[2] # 導波路損失(dB/cm)
gamc = np.log(10**(gac*1E-5)) # 導波路損失(/um)
Sc = [0 for a in range(div)]
Sfc = [0 for a in range(div)]
Smaxc = 0

i = 0
while i < div:
    Sc[i] = 1 - np.exp(-gamc*z[i]) + np.exp(-gamc*z[i])*np.exp(-2*g*z[i])#(gamc + 2*g*np.exp(-(gamc+2*g)*z[i]))/(gamc+2*g)
    Sfc[i] = -10*np.log10(Sc[i])
    if Smaxc < Sfc[i]:
        Smaxc = Sfc[i]
        zmaxc = z[i]
    i += 1

#-----------------------------------------------------------------------------
# 3
gad = a[3] # 導波路損失(dB/cm)
gamd = np.log(10**(gad*1E-5)) # 導波路損失(/um)
Sd = [0 for a in range(div)]
Sfd = [0 for a in range(div)]
Smaxd = 0

i = 0
while i < div:
    Sd[i] = 1 - np.exp(-gamd*z[i]) + np.exp(-gamd*z[i])*np.exp(-2*g*z[i])#(gamd + 2*g*np.exp(-(gamd+2*g)*z[i]))/(gamd+2*g)
    Sfd[i] = -10*np.log10(Sd[i])
    if Smaxd < Sfd[i]:
        Smaxd = Sfd[i]
        zmaxd = z[i]
    i += 1

#-----------------------------------------------------------------------------
# 4
gae = a[4] # 導波路損失(dB/cm)
game = np.log(10**(gae*1E-5)) # 導波路損失(/um)
Se = [0 for a in range(div)]
Sfe = [0 for a in range(div)]
Smaxe = 0

i = 0
while i < div:
    Se[i] =  1 - np.exp(-game*z[i]) + np.exp(-game*z[i])*np.exp(-2*g*z[i]) #(game + 2*g*np.exp(-(game+2*g)*z[i]))/(game+2*g)
    Sfe[i] = -10*np.log10(Se[i])
    if Smaxe < Sfe[i]:
        Smaxe = Sfe[i]
        zmaxe = z[i]
    i += 1


print("Smax = {} dB (loss = {} dB/cm) Zmax = {} um".format(Smaxa, a[0], zmaxa))
print("Smax = {} dB (loss = {} dB/cm) Zmax = {} um".format(Smaxb, a[1], zmaxb))
print("Smax = {} dB (loss = {} dB/cm) Zmax = {} um".format(Smaxc, a[2], zmaxc))
print("Smax = {} dB (loss = {} dB/cm) Zmax = {} um".format(Smaxd, a[3], zmaxd))
print("Smax = {} dB (loss = {} dB/cm) Zmax = {} um".format(Smaxe, a[4], zmaxe))
#グラフを描く
plt.plot(z/1000, Sfa, color='purple', label=str(a[0])+' dB/cm')
plt.plot(z/1000, Sfb, color='red', label=str(a[1])+' dB/cm')
plt.plot(z/1000, Sfc, color='orange', label=str(a[2])+' dB/cm')
plt.plot(z/1000, Sfd, color='green', label=str(a[3])+' dB/cm')
plt.plot(z/1000, Sfe, color='blue', label=str(a[4])+' dB/cm')
plt.xlabel('Device length (mm)')
plt.ylabel('Squeezing level (dB)')
plt.legend(loc='upper left', title='Propagation loss')#loc='best')
plt.gca().xaxis.set_tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)
plt.gca().yaxis.set_tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)
plt.xlim([0, 40])
plt.ylim([0, 10])

plt.show()
