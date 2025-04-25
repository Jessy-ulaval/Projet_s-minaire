import sys

sys.path.append("/fonctions_projet.py")

from fonctions_projet import *
import matplotlib.pyplot as plt


############################# Situation 2 ##################################

# Tarifs optimaux
s_t_o_b_2_vect, s_t_o_b_2_dict = tarifs_optimaux_bilateraux(1,1, 1, c, c, c, c, c, c)


s_t_o_m_2_vect, s_t_o_m_2_dict = tarifs_optimaux_multilateraux(1,1, 1, c, c, c, c, c, c)


# Prix d'équilibre
p_equ_vect, p_equ_dict = prix_equilibre(a, 1, 1, 1, s_t_o_b_2_dict["ti_quo"], s_t_o_b_2_dict["tj_quo"],
                                        s_t_o_b_2_dict["tk_quo"], s_t_o_b_2_dict["ti_fij"], s_t_o_b_2_dict["tj_fij"]
                                        , s_t_o_b_2_dict["ti_fik"], s_t_o_b_2_dict["tk_fik"], s_t_o_b_2_dict["tj_fjk"],
                                        s_t_o_b_2_dict["tk_fjk"], s_t_o_m_2_dict["ti_mij"]
                                        , s_t_o_m_2_dict["tj_mij"], s_t_o_m_2_dict["ti_mik"], s_t_o_m_2_dict["tk_mik"],
                                        s_t_o_m_2_dict["tj_mjk"], s_t_o_m_2_dict["tk_mjk"]
                                        , s_t_o_b_2_dict["tj_ih"], s_t_o_b_2_dict["tk_ih"], s_t_o_b_2_dict["ti_jh"],
                                        s_t_o_b_2_dict["tk_jh"], s_t_o_b_2_dict["ti_kh"], s_t_o_b_2_dict["tj_kh"]
                                        , c, c, c, c, c, c)

# Exportation
x_quo = export_i(a, 1, 1, 1, p_equ_dict["pi_i_quo"], p_equ_dict["pj_j_quo"], p_equ_dict["pk_k_quo"],
                 s_t_o_b_2_dict["ti_quo"], s_t_o_b_2_dict["ti_quo"], s_t_o_b_2_dict["tj_quo"], s_t_o_b_2_dict["tj_quo"],
                 s_t_o_b_2_dict["tk_quo"], s_t_o_b_2_dict["tk_quo"], c, c, c, c, c, c)
x_fij = export_i(a, 1, 1, 1, p_equ_dict["pi_i_fij"], p_equ_dict["pj_j_fij"], p_equ_dict["pk_k_fij"], 0,
                 s_t_o_b_2_dict["ti_fij"], 0, s_t_o_b_2_dict["tj_fij"], s_t_o_b_2_dict["tk_quo"],
                 s_t_o_b_2_dict["tk_quo"], c, c, c, c, c, c)
x_fik = export_i(a, 1, 1, 1, p_equ_dict["pi_i_fik"], p_equ_dict["pj_j_fik"], p_equ_dict["pk_k_fik"],
                 s_t_o_b_2_dict["ti_fik"], 0, s_t_o_b_2_dict["tj_quo"], s_t_o_b_2_dict["tj_quo"], 0,
                 s_t_o_b_2_dict["tk_fik"], c, c, c, c, c, c)
x_fjk = export_i(a, 1, 1, 1, p_equ_dict["pi_i_fjk"], p_equ_dict["pj_j_fjk"], p_equ_dict["pk_k_fjk"],
                 s_t_o_b_2_dict["ti_quo"], s_t_o_b_2_dict["ti_quo"], s_t_o_b_2_dict["tj_fjk"], 0,
                 s_t_o_b_2_dict["tk_fjk"], 0, c, c, c, c, c, c)
x_mij = export_i(a, 1, 1, 1, p_equ_dict["pi_i_mij"], p_equ_dict["pj_j_mij"], p_equ_dict["pk_k_mij"],
                 s_t_o_m_2_dict["ti_mij"], s_t_o_m_2_dict["ti_mij"], s_t_o_m_2_dict["tj_mij"], s_t_o_m_2_dict["tj_mij"],
                 s_t_o_m_2_dict["tk_quo"], s_t_o_m_2_dict["tk_quo"], c, c, c, c, c, c)
x_mik = export_i(a, 1, 1, 1, p_equ_dict["pi_i_mik"], p_equ_dict["pj_j_mik"], p_equ_dict["pk_k_mik"],
                 s_t_o_m_2_dict["ti_mik"], s_t_o_m_2_dict["ti_mik"], s_t_o_m_2_dict["tj_quo"], s_t_o_m_2_dict["tj_quo"],
                 s_t_o_m_2_dict["tk_mik"], s_t_o_m_2_dict["tk_mik"], c, c, c, c, c, c)
x_mjk = export_i(a, 1, 1, 1, p_equ_dict["pi_i_mjk"], p_equ_dict["pj_j_mjk"], p_equ_dict["pk_k_mjk"],
                 s_t_o_m_2_dict["ti_quo"], s_t_o_m_2_dict["ti_quo"], s_t_o_m_2_dict["tj_mjk"], s_t_o_m_2_dict["tj_mjk"],
                 s_t_o_m_2_dict["tk_mjk"], s_t_o_m_2_dict["tk_mjk"], c, c, c, c, c, c)
x_ih = export_i(a, 1, 1, 1, p_equ_dict["pi_i_ih"], p_equ_dict["pj_j_ih"], p_equ_dict["pk_k_ih"], 0, 0,
                0,
                s_t_o_b_2_dict["tj_ih"], 0, s_t_o_b_2_dict["tk_ih"], c, c, c, c, c, c)
x_jh = export_i(a, 1, 1, 1, p_equ_dict["pi_i_jh"], p_equ_dict["pj_j_jh"], p_equ_dict["pk_k_jh"], 0,
                s_t_o_b_2_dict["ti_jh"], 0, 0, s_t_o_b_2_dict["tk_jh"], 0, c, c, c, c, c, c)
x_kh = export_i(a, 1, 1, 1, p_equ_dict["pi_i_kh"], p_equ_dict["pj_j_kh"], p_equ_dict["pk_k_kh"],
                s_t_o_b_2_dict["ti_kh"], 0, s_t_o_b_2_dict["tj_kh"], 0, 0, 0, c, c, c, c, c, c)
x_F = export_i(a, 1, 1, 1, p_equ_dict["pi_i_F"], p_equ_dict["pj_j_F"], p_equ_dict["pk_k_F"], 0, 0,
               0, 0, 0, 0, c, c, c, c, c, c)
X_b, X_m = export_final(x_quo, x_fij, x_fik, x_fjk, x_mij, x_mik, x_mjk, x_ih, x_jh, x_kh, x_F)


# Borne c
obtenir_bornes_restrictives(s_t_o_b_2_vect, s_t_o_b_2_dict, s_t_o_m_2_vect, s_t_o_m_2_dict, X_b, X_m, c)
print("")


# Bien-être optimaux - Structure bilaterale
print("##### Statu quo ####")
w_quo = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, s_t_o_b_2_dict["ti_quo"], s_t_o_b_2_dict["tj_quo"],
                   s_t_o_b_2_dict["ti_quo"], s_t_o_b_2_dict["tk_quo"], s_t_o_b_2_dict["tj_quo"],
                   s_t_o_b_2_dict["tk_quo"],"quo")
print("")

print("#### fij ####")
w_fij = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, 0, 0, s_t_o_b_2_dict["ti_fij"],
                   s_t_o_b_2_dict["tk_quo"],
                   s_t_o_b_2_dict["tj_fij"], s_t_o_b_2_dict["tk_quo"], "fij")
print("")

print("##### fik ####")
w_fik = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, s_t_o_b_2_dict["ti_fik"], s_t_o_b_2_dict["tj_quo"], 0,
                   0,
                   s_t_o_b_2_dict["tj_quo"], s_t_o_b_2_dict["tk_fik"], "fik")
print("")

print("##### fjk ####")
w_fjk = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, s_t_o_b_2_dict["ti_quo"], s_t_o_b_2_dict["tj_fjk"],
                   s_t_o_b_2_dict["ti_quo"], s_t_o_b_2_dict["tk_fjk"], 0, 0, "fjk")
print("")

print("##### ih ####")
w_ih = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, 0, 0, 0, 0, s_t_o_b_2_dict["tj_ih"],
                  s_t_o_b_2_dict["tk_ih"], "ih")
print("")

print("##### jh ####")
w_jh = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, 0, 0, s_t_o_b_2_dict["ti_jh"], s_t_o_b_2_dict["tk_jh"],
                  0, 0, "jh")
print("")

print("##### kh ####")
w_kh = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, s_t_o_b_2_dict["ti_kh"], s_t_o_b_2_dict["tj_kh"], 0,
                  0, 0, 0, "kh")
print("")

print("##### F ####")
w_F = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, 0, 0, 0, 0, 0, 0, "F")
print("")

# Bien-être optimaux - Structure multilateral
print("##### mij ####")
w_mij = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, s_t_o_m_2_dict["ti_mij"], s_t_o_m_2_dict["tj_mij"],
                   s_t_o_m_2_dict["ti_mij"], s_t_o_m_2_dict["tk_quo"], s_t_o_m_2_dict["tj_mij"],
                   s_t_o_m_2_dict["tk_quo"], "mij")
print("")

print("##### mik ####")
w_mik = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, s_t_o_m_2_dict["ti_mik"], s_t_o_m_2_dict["tj_quo"],
                   s_t_o_m_2_dict["ti_mik"], s_t_o_m_2_dict["tk_mik"], s_t_o_m_2_dict["tj_quo"],
                   s_t_o_m_2_dict["tk_mik"], "mik")
print("")

print("##### mjk ####")
w_mjk = w_optimaux(a, 1, 1, 1, c, c, c, c, c, c, s_t_o_m_2_dict["ti_quo"], s_t_o_m_2_dict["tj_mjk"],
                   s_t_o_m_2_dict["ti_quo"], s_t_o_m_2_dict["tk_mjk"], s_t_o_m_2_dict["tj_mjk"],
                   s_t_o_m_2_dict["tk_mjk"], "mjk")
print("")





######################################### Déviation unilatérale - Bilatérale #########################################
c_values_b = np.arange(0, 1, 0.001)
n_b = len(c_values_b)
vect_pays = ['i','j','k']
structure = "bilatérale sys"

# Déviation unilatérale - quo
EN_quo = []
for i in range(n_b):
    EN_quo.append(1)

# Déviation unilatérale - fij
vect_dev_fij = [w_quo]
symbole_dev = ['quo']
EN_fij = []
deviation_fij, dev_inverse_fij, dev_global_fij = deviation_uni("fij",symbole_dev, w_fij, vect_dev_fij,
                                                               vect_pays, structure, n_b, c_values_b)
for i in range(n_b):
    if deviation_fij['wi_fij_quo'][i] == 1 and deviation_fij['wj_fij_quo'][i] == 1:
        EN_fij.append(1)
    else: EN_fij.append(0)
print(EN_fij)
inverse_fij = [1 - x for x in EN_fij]

fig, ax = plt.subplots()
axe(structure,deviation_fij["wi_fij_quo"], ax, c_values_b, 0)
axe(structure,deviation_fij["wj_fij_quo"], ax, c_values_b, 1)
axe(structure,EN_fij, ax, c_values_b, 2)
ax.set_yticks([0, 1, 2])
ax.set_yticklabels(["$w_i$($f_{ij}$-$\phi$)", "$w_j$($f_{ij}$-$\phi$)", "EN_$f_{ij}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale bilatérale $f_{ij}$")
plt.show()



# Déviation unilatérale - fik
vect_dev_fik = [w_quo]
symbole_dev = ['quo']
EN_fik = []
deviation_fik, dev_inverse_fik, dev_global_fik = deviation_uni("fik",symbole_dev, w_fik, vect_dev_fik,
                                                               vect_pays, structure, n_b, c_values_b)
for i in range(n_b):
    if deviation_fik['wi_fik_quo'][i] == 1 and deviation_fik['wk_fik_quo'][i] == 1:
        EN_fik.append(1)
    else: EN_fik.append(0)
print(EN_fik)
inverse_fik = [1 - x for x in EN_fik]

fig, ax = plt.subplots()
axe(structure,deviation_fik["wi_fik_quo"], ax, c_values_b, 0)
axe(structure,deviation_fik["wk_fik_quo"], ax, c_values_b, 1)
axe(structure,EN_fik, ax, c_values_b, 2)
ax.set_yticks([0, 1, 2])
ax.set_yticklabels(["$w_i$($f_{ik}$-$\phi$)", "$w_k$($f_{ik}$-$\phi$)", "EN_$f_{ik}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale bilatérale $f_{ik}$")
plt.show()

# Déviation unilatérale - fjk
vect_dev_fjk = [w_quo]
symbole_dev = ['quo']
EN_fjk = []
deviation_fjk, dev_inverse_fjk, dev_global_fjk = deviation_uni("fjk",symbole_dev, w_fjk, vect_dev_fjk,
                                                               vect_pays, structure, n_b, c_values_b)
for i in range(n_b):
    if deviation_fjk['wj_fjk_quo'][i] == 1 and deviation_fjk['wk_fjk_quo'][i] == 1:
        EN_fjk.append(1)
    else: EN_fjk.append(0)
print(EN_fjk)
inverse_fjk = [1 - x for x in EN_fjk]

fig, ax = plt.subplots()
axe(structure,deviation_fjk["wj_fjk_quo"], ax, c_values_b, 0)
axe(structure,deviation_fjk["wk_fjk_quo"], ax, c_values_b, 1)
axe(structure,EN_fjk, ax, c_values_b, 2)
ax.set_yticks([0, 1, 2])
ax.set_yticklabels(["$w_j$($f_{jk}$-$\phi$)", "$w_k$($f_{jk}$-$\phi$)", "EN_$f_{jk}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale bilatérale $f_{jk}$")
plt.show()







# Déviation unilatérale - ih
vect_dev_ih = [w_fij, w_fik, w_quo]
symbole_dev = ['fij', 'fik', 'quo']
EN_ih = []
deviation_ih, dev_inverse_ih, dev_global_ih = deviation_uni("ih",symbole_dev, w_ih, vect_dev_ih,
                                                            vect_pays, structure, n_b, c_values_b)
for i in range(n_b):
    if (deviation_ih['wi_ih_fij'][i] == 1 and deviation_ih['wi_ih_fik'][i] == 1 and deviation_ih['wi_ih_quo'][i] == 1
            and deviation_ih['wj_ih_fik'][i] == 1 and deviation_ih['wk_ih_fij'][i] == 1):
        EN_ih.append(1)
    else: EN_ih.append(0)
print(EN_ih)
inverse_ih = [1 - x for x in EN_ih]

fig, ax = plt.subplots()
axe(structure,deviation_ih["wi_ih_fij"],ax,c_values_b,0)
axe(structure,deviation_ih["wi_ih_fik"], ax, c_values_b, 1)
axe(structure,deviation_ih["wi_ih_quo"], ax, c_values_b, 2)
axe(structure,deviation_ih["wj_ih_fik"], ax, c_values_b, 3)
axe(structure,deviation_ih["wk_ih_fij"], ax, c_values_b, 4)
axe(structure,EN_ih, ax, c_values_b, 5)
ax.set_yticks([0, 1, 2, 3, 4, 5])
ax.set_yticklabels(["$w_i$(ih-$f_{ij}$)","$w_i$(ih-$f_{ik}$)","$w_i$(ih-$\phi$)","$w_j$(ih-$f_{ik}$)",
                    "$w_k$(ih-$f_{ij}$)", "EN_ih"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale bilatérale ih")
plt.show()

# Déviation unilatérale - jh
vect_dev_jh = [w_fjk, w_fij, w_quo]
symbole_dev = ['fjk', 'fij', 'quo']
EN_jh = []
deviation_jh, dev_inverse_jh, dev_global_jh= deviation_uni("jh",symbole_dev, w_jh, vect_dev_jh, vect_pays,
                                                           structure, n_b, c_values_b)
for i in range(n_b):
    if (deviation_jh['wj_jh_fjk'][i] == 1 and deviation_jh['wj_jh_fij'][i] == 1 and deviation_jh['wj_jh_quo'][i] == 1
            and deviation_jh['wk_jh_fij'][i] == 1 and deviation_jh['wi_jh_fjk'][i] == 1):
        EN_jh.append(1)
    else: EN_jh.append(0)
print(EN_jh)
inverse_jh = [1 - x for x in EN_jh]

fig, ax = plt.subplots()
axe(structure,deviation_jh["wj_jh_fjk"],ax,c_values_b,0)
axe(structure,deviation_jh["wj_jh_fij"], ax, c_values_b, 1)
axe(structure,deviation_jh["wj_jh_quo"], ax, c_values_b, 2)
axe(structure,deviation_jh["wk_jh_fij"], ax, c_values_b, 3)
axe(structure,deviation_jh["wi_jh_fjk"], ax, c_values_b, 4)
axe(structure,EN_jh, ax, c_values_b, 5)
ax.set_yticks([0, 1, 2, 3, 4, 5])
ax.set_yticklabels(["$w_j$(jh-$f_{jk}$)","$w_j$(jh-$f_{ij}$)","$w_j$(jh-$\phi$)","$w_k$(jh-$f_{ij}$)",
                    "$w_i$(jh-$f_{jk}$)", "EN_jh"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale bilatérale jh")
plt.show()



# Déviation unilatérale - kh
vect_dev_kh = [w_fjk, w_fik, w_quo]
symbole_dev = ['fjk', 'fik', 'quo']
EN_kh = []
deviation_kh, dev_inverse_kh, dev_global_kh = deviation_uni("kh",symbole_dev, w_kh, vect_dev_kh, vect_pays,
                                                            structure, n_b, c_values_b)
for i in range(n_b):
    if (deviation_kh['wk_kh_fjk'][i] == 1 and deviation_kh['wk_kh_fik'][i] == 1 and deviation_kh['wk_kh_quo'][i] == 1
            and deviation_kh['wj_kh_fik'][i] == 1 and deviation_kh['wi_kh_fjk'][i] == 1):
        EN_kh.append(1)
    else: EN_kh.append(0)
print(EN_kh)
inverse_kh = [1 - x for x in EN_kh]

fig, ax = plt.subplots()
axe(structure,deviation_kh["wk_kh_fjk"],ax,c_values_b,0)
axe(structure,deviation_kh["wk_kh_fik"], ax, c_values_b, 1)
axe(structure,deviation_kh["wk_kh_quo"], ax, c_values_b, 2)
axe(structure,deviation_kh["wj_kh_fik"], ax, c_values_b, 3)
axe(structure,deviation_kh["wi_kh_fjk"], ax, c_values_b, 4)
axe(structure,EN_kh, ax, c_values_b, 5)
ax.set_yticks([0, 1, 2, 3, 4, 5])
ax.set_yticklabels(["$w_k$(kh-$f_{jk}$)","$w_k$(kh-$f_{ik}$)","$w_k$(kh-$\phi$)","$w_j$(kh-$f_{ik}$)",
                    "$w_i$(kh-$f_{jk}$)", "EN_kh"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale bilatérale kh")
plt.show()

# Déviation unilatérale - F
vect_dev_F = [w_jh, w_kh, w_fjk, w_ih, w_fik, w_fij]
symbole_dev = ['jh', 'kh', 'fjk', 'ih', 'fik', 'fij']
EN_F = []
deviation_F, dev_inverse_F, dev_global_F = deviation_uni("F",symbole_dev, w_F, vect_dev_F, vect_pays,
                                                         structure, n_b, c_values_b)
for i in range(n_b):
    if (deviation_F['wi_F_jh'][i] == 1 and deviation_F['wi_F_kh'][i] == 1 and deviation_F['wi_F_fjk'][i] == 1 and
            deviation_F['wj_F_ih'][i] == 1 and deviation_F['wj_F_kh'][i] == 1 and deviation_F['wj_F_fik'][i] == 1 and
                deviation_F['wk_F_ih'][i] == 1 and deviation_F['wk_F_jh'][i] == 1 and deviation_F['wk_F_fij'][i] == 1):
        EN_F.append(1)
    else: EN_F.append(0)
print(EN_F)
inverse_F = [1 - x for x in EN_F]

fig, ax = plt.subplots()
axe(structure,deviation_F["wi_F_jh"],ax,c_values_b,0)
axe(structure,deviation_F["wi_F_kh"],ax,c_values_b,1)
axe(structure,deviation_F["wi_F_fjk"],ax,c_values_b,2)
axe(structure,deviation_F["wj_F_ih"],ax,c_values_b,3)
axe(structure,deviation_F["wj_F_kh"],ax,c_values_b,4)
axe(structure,deviation_F["wj_F_fik"],ax,c_values_b,5)
axe(structure,deviation_F["wk_F_ih"],ax,c_values_b,6)
axe(structure,deviation_F["wk_F_jh"],ax,c_values_b,7)
axe(structure,deviation_F["wk_F_fij"],ax,c_values_b,8)
axe(structure,EN_F,ax,c_values_b,9)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
ax.set_yticklabels(["$w_i$(F-jh)","$w_i$(F-kh)","$w_i$(F-$f_{jk}$)","$w_j$(F-ih)", "$w_j$(F-kh)", "$w_j$(F-$f_{ik}$)",
                    "$w_k$(F-ih)", "$w_k$(F-jh)", "$w_k$(F-$f_{ij}$)", "EN_F"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale bilatérale F")
plt.show()

# Graphique EN bilatérale
fig, ax = plt.subplots()
axe(structure,EN_quo, ax, c_values_b, 0)
axe(structure,EN_fij, ax, c_values_b, 1)
axe(structure,EN_fik, ax, c_values_b, 2)
axe(structure,EN_fjk, ax, c_values_b, 3)
axe(structure,EN_ih, ax, c_values_b, 4)
axe(structure,EN_jh, ax, c_values_b, 5)
axe(structure,EN_kh, ax, c_values_b, 6)
axe(structure,EN_F, ax, c_values_b, 7)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_yticklabels(["EN_$\phi$", "EN_$f_{ij}$", "EN_$f_{ik}$", "EN_$f_{jk}$", "EN_ih", "EN_jh", "EN_kh", "EN_F"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Équilibre de Nash - Déviation unilatérale bilatérale")
plt.show()


######################################### Déviation unilatérale - multilatérale ###################################
c_values_m = np.arange(0, 1, 0.001)
n_m = len(c_values_m)
structure = "multilatérale"
# Déviation unilatérale - quo
EN_m_quo = []
for i in range(n_m):
    EN_m_quo.append(1)

# Déviation unilatérale - mij
vect_dev_mij = [w_quo]
symbole_dev = ['quo']
EN_mij = []
deviation_mij, dev_inverse_mij, dev_global_mij = deviation_uni("mij",symbole_dev, w_mij, vect_dev_mij,
                                                               vect_pays, structure, n_m, c_values_m)
for i in range(n_m):
    if deviation_mij['wi_mij_quo'][i] == 1 and deviation_mij['wj_mij_quo'][i] == 1:
        EN_mij.append(1)
    else: EN_mij.append(0)
print(EN_mij)
inverse_mij = [1 - x for x in EN_mij]

fig, ax = plt.subplots()
axe(structure,deviation_mij["wi_mij_quo"], ax, c_values_m, 0)
axe(structure,deviation_mij["wj_mij_quo"], ax, c_values_m, 1)
axe(structure,EN_mij, ax, c_values_m, 2)
ax.set_yticks([0, 1, 2])
ax.set_yticklabels(["$w_i$($m_{ij}$-$\phi$)", "$w_j$($m_{ij}$-$\phi$)", "EN_$m_{ij}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale multilatérale $m_{ij}$")
plt.show()

# Déviation unilatérale - mik
vect_dev_mik = [w_quo]
symbole_dev = ['quo']
EN_mik = []
deviation_mik, dev_inverse_mik, dev_global_mik = deviation_uni("mik",symbole_dev, w_mik, vect_dev_mik,
                                                               vect_pays, structure, n_m, c_values_m)
for i in range(n_m):
    if deviation_mik['wi_mik_quo'][i] == 1 and deviation_mik['wk_mik_quo'][i] == 1:
        EN_mik.append(1)
    else: EN_mik.append(0)
print(EN_mik)
inverse_mik = [1 - x for x in EN_mik]

fig, ax = plt.subplots()
axe(structure,deviation_mik["wi_mik_quo"], ax, c_values_m, 0)
axe(structure,deviation_mik["wk_mik_quo"], ax, c_values_m, 1)
axe(structure,EN_mik, ax, c_values_m, 2)
ax.set_yticks([0, 1, 2])
ax.set_yticklabels(["$w_i$($m_{ik}$-$\phi$)", "$w_k$($m_{ik}$-$\phi$)", "EN_$m_{ik}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale multilatérale $m_{ik}$")
plt.show()

# Déviation unilatérale - mjk
vect_dev_mjk = [w_quo]
symbole_dev = ['quo']
EN_mjk = []
deviation_mjk, dev_inverse_mjk, dev_global_mjk = deviation_uni("mjk",symbole_dev, w_mjk, vect_dev_mjk,
                                                               vect_pays, structure, n_m, c_values_m)
for i in range(n_m):
    if deviation_mjk['wj_mjk_quo'][i] == 1 and deviation_mjk['wk_mjk_quo'][i] == 1:
        EN_mjk.append(1)
    else: EN_mjk.append(0)
print(EN_mjk)
inverse_mjk = [1 - x for x in EN_mjk]

fig, ax = plt.subplots()
axe(structure,deviation_mjk["wj_mjk_quo"], ax, c_values_m, 0)
axe(structure,deviation_mjk["wk_mjk_quo"], ax, c_values_m, 1)
axe(structure,EN_mjk, ax, c_values_m, 2)
ax.set_yticks([0, 1, 2])
ax.set_yticklabels(["$w_j$($m_{jk}$-$\phi$)", "$w_k$($m_{jk}$-$\phi$)", "EN_$m_{jk}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unlatérale multilatérale $m_{jk}$")
plt.show()


# Déviation unilatérale - F
vect_dev_F = [w_mjk, w_mik, w_mij]
symbole_dev = ['mjk', 'mik', 'mij']
EN_m_F = []
deviation_m_F, dev_inverse_m_F, dev_global_m_F = deviation_uni("F",symbole_dev, w_F, vect_dev_F, vect_pays,
                                                               structure, n_m, c_values_m)
for i in range(n_m):
    if (deviation_m_F['wi_F_mjk'][i] == 1 and deviation_m_F['wj_F_mik'][i] == 1 and deviation_m_F['wk_F_mij'][i] == 1):
        EN_m_F.append(1)
    else: EN_m_F.append(0)
print(EN_m_F)
inverse_m_F = [1 - x for x in EN_m_F]

ig, ax = plt.subplots()
axe(structure,deviation_m_F["wi_F_mjk"],ax,c_values_m,0)
axe(structure,deviation_m_F["wj_F_mik"],ax,c_values_m,1)
axe(structure,deviation_m_F["wk_F_mij"],ax,c_values_m,2)
axe(structure,EN_m_F,ax,c_values_m,3)
ax.set_yticks([0, 1, 2, 3])
ax.set_yticklabels(["$w_i$(F-$m_{jk}$)","$w_j$(F-$m_{jk}$)", "$w_k$(F-$m_{ij}$)", "EN_F"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation unilatérale multilatérale F")
plt.show()


# Graphique EN multilatérale
fig, ax = plt.subplots()
axe(structure,EN_m_quo, ax, c_values_m, 0)
axe(structure,EN_mij, ax, c_values_m, 1)
axe(structure,EN_mik, ax, c_values_m, 2)
axe(structure,EN_mjk, ax, c_values_m, 3)
axe(structure,EN_m_F, ax, c_values_m, 4)
ax.set_yticks([0, 1, 2, 3, 4])
ax.set_yticklabels(["EN_$\phi$", "EN_$m_{ij}$", "EN_$m_{ik}$", "EN_$m_{jk}$", "EN_F"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Équilibre de Nash - Déviation unilatérale multilatérale")

plt.show()


##################################### Déviation de coalition - Bilatérale ##############################################
structure = "bilatérale sys"
# Déviation coalition - F
vect_dev_F = [w_quo, w_quo, w_quo, w_fij, w_fik, w_fjk]
symbole_dev = ['quo', 'quo', 'quo', 'fij', 'fik', 'fjk']
vect_pays = [['i','j'], ['i','k'], ['j','k'],['i', 'j'],['i','k'], ['j','k']]
vect_cre_F = [[[dev_global_fik['w_fik_quo'][0]],[dev_global_fjk['w_fjk_quo'][1]]],
              [[dev_global_fij['w_fij_quo'][0]],[dev_global_fjk['w_fjk_quo'][2]]],
              [[dev_global_fij['w_fij_quo'][1]],[dev_global_fik['w_fik_quo'][2]]],
              [[dev_global_ih['w_ih_fij'][0],dev_inverse_fij['w_quo_fij'][0]], [dev_global_jh['w_jh_fij'][1],
                                                                                dev_inverse_fij['w_quo_fij'][1]]],
              [[dev_global_ih['w_ih_fik'][0],dev_inverse_fik['w_quo_fik'][0]], [dev_global_kh['w_kh_fik'][2],
                                                                                dev_inverse_fik['w_quo_fik'][2]]],
              [[dev_global_jh['w_jh_fjk'][1],dev_inverse_fjk['w_quo_fjk'][1]], [dev_global_kh['w_kh_fjk'][2],
                                                                                dev_inverse_fjk['w_quo_fjk'][2]]],
              ]
print("")


ENEC_F = []
deviation_enec_F = deviation_coal("F",symbole_dev, w_F, vect_dev_F, vect_cre_F, vect_pays, structure, n_b,
                                  c_values_b)
for i in range(n_b):
    if EN_F[i] == 0:
        ENEC_F.append(0)
    elif ((deviation_enec_F['wi_F_quo_1'][i] + deviation_enec_F['wj_F_quo_1'][i]) == 0
          or (deviation_enec_F['wi_F_quo_2'][i] + deviation_enec_F['wk_F_quo_2'][i]) == 0
          or (deviation_enec_F['wj_F_quo_3'][i] + deviation_enec_F['wk_F_quo_3'][i]) == 0
          or (deviation_enec_F['wi_F_fij_4'][i] + deviation_enec_F['wj_F_fij_4'][i]) == 0
          or (deviation_enec_F['wi_F_fik_5'][i] + deviation_enec_F['wk_F_fik_5'][i]) == 0
          or (deviation_enec_F['wj_F_fjk_6'][i] + deviation_enec_F['wk_F_fjk_6'][i]) == 0):
        ENEC_F.append(0)
    else: ENEC_F.append(1)
print(ENEC_F)
inverse_ENEC_F = [1 - x for x in ENEC_F]

fig, ax = plt.subplots()
axe(structure,deviation_enec_F["wi_F_quo_1"], ax, c_values_b, 0)
axe(structure,deviation_enec_F["wj_F_quo_1"], ax, c_values_b, 1)
axe(structure,deviation_enec_F["wi_F_quo_2"], ax, c_values_b, 2)
axe(structure,deviation_enec_F["wk_F_quo_2"], ax, c_values_b, 3)
axe(structure,deviation_enec_F["wj_F_quo_3"], ax, c_values_b, 4)
axe(structure,deviation_enec_F["wk_F_quo_3"], ax, c_values_b, 5)
axe(structure,deviation_enec_F["wi_F_fij_4"], ax, c_values_b, 6)
axe(structure,deviation_enec_F["wj_F_fij_4"], ax, c_values_b, 7)
axe(structure,deviation_enec_F["wi_F_fik_5"], ax, c_values_b, 8)
axe(structure,deviation_enec_F["wk_F_fik_5"], ax, c_values_b, 9)
axe(structure,deviation_enec_F["wj_F_fjk_6"], ax, c_values_b, 10)
axe(structure,deviation_enec_F["wk_F_fjk_6"], ax, c_values_b, 11)
axe(structure,ENEC_F, ax, c_values_b, 12)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6,7,8,9,10,11,12])
ax.set_yticklabels(["$w_i$(F-$\phi$)", "$w_j$(F-$\phi$)","$w_i$(F-$\phi$)", "$w_k$(F-$\phi$)", "$w_j$(F-$\phi$)",
                    "$w_k$(F-$\phi$)",
                    "$w_i$(F-$f_{ij}$)", "$w_j$(F-$f_{ij}$)", "$w_i$(F-$f_{ik}$)","$w_k$(F-$f_{ik}$)",
                    "$w_j$(F-$f_{jk}$)", "$w_k$(F-$f_{jk}$)", "ENEC_F"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition bilatérale F")
plt.show()




# Déviation coalition - statu quo
vect_dev_quo = [w_fij, w_fik, w_fjk, w_ih, w_jh, w_kh, w_F]
symbole_dev = ['fij', 'fik', 'fjk', 'ih', 'jh', 'kh', 'F']
vect_pays = [['i','j'], ['i','k'], ['j','k'],['i', 'j', 'k'],['i', 'j', 'k'], ['i', 'j', 'k'], ['i', 'j', 'k']]
vect_cre_quo = [[[dev_inverse_fij['w_quo_fij'][0]],[dev_inverse_fij['w_quo_fij'][1]]],
                [[dev_inverse_fik['w_quo_fik'][0]],[dev_inverse_fik['w_quo_fik'][2]]],
                [[dev_inverse_fjk['w_quo_fjk'][1]],[dev_inverse_fjk['w_quo_fjk'][2]]],
                [[inverse_ih], [inverse_ih, dev_global_F['w_F_ih'][1]], [inverse_ih, dev_global_F['w_F_ih'][2]]],
                [[inverse_jh, dev_global_F['w_F_jh'][0]], [inverse_jh],  [inverse_jh, dev_global_F['w_F_jh'][2]]],
                [[inverse_kh, dev_global_F['w_F_kh'][0]], [inverse_kh, dev_global_F['w_F_kh'][1]], [inverse_kh]],
                [[inverse_ENEC_F], [inverse_ENEC_F],[inverse_ENEC_F]]
                ]


ENEC_quo = []
deviation_enec_quo = deviation_coal('quo', symbole_dev, w_quo, vect_dev_quo, vect_cre_quo, vect_pays,
                                    structure, n_b, c_values_b)
for i in range(n_b):
    if EN_quo[i] == 0:
        ENEC_quo.append(0)
    elif ((deviation_enec_quo['wi_quo_fij_1'][i] + deviation_enec_quo['wj_quo_fij_1'][i]) == 0
        or (deviation_enec_quo['wi_quo_fik_2'][i] + deviation_enec_quo['wk_quo_fik_2'][i]) == 0
        or (deviation_enec_quo['wj_quo_fjk_3'][i] + deviation_enec_quo['wk_quo_fjk_3'][i]) == 0
        or (deviation_enec_quo['wi_quo_ih_4'][i] + deviation_enec_quo['wj_quo_ih_4'][i] +
            deviation_enec_quo['wk_quo_ih_4'][i]) == 0
        or (deviation_enec_quo['wi_quo_jh_5'][i] + deviation_enec_quo['wj_quo_jh_5'][i] +
            deviation_enec_quo['wk_quo_jh_5'][i]) == 0
        or (deviation_enec_quo['wi_quo_kh_6'][i] + deviation_enec_quo['wj_quo_kh_6'][i] +
            deviation_enec_quo['wk_quo_kh_6'][i]) == 0
        or (deviation_enec_quo['wi_quo_F_7'][i] + deviation_enec_quo['wj_quo_F_7'][i] +
            deviation_enec_quo['wk_quo_F_7'][i]) == 0):
        ENEC_quo.append(0)
    else: ENEC_quo.append(1)
print(ENEC_quo)

fig, ax = plt.subplots()
axe(structure,deviation_enec_quo["wi_quo_fij_1"], ax, c_values_b, 0)
axe(structure,deviation_enec_quo["wj_quo_fij_1"], ax, c_values_b, 1)
axe(structure,deviation_enec_quo["wi_quo_fik_2"], ax, c_values_b, 2)
axe(structure,deviation_enec_quo["wk_quo_fik_2"], ax, c_values_b, 3)
axe(structure,deviation_enec_quo["wj_quo_fjk_3"], ax, c_values_b, 4)
axe(structure,deviation_enec_quo["wk_quo_fjk_3"], ax, c_values_b, 5)
axe(structure,deviation_enec_quo["wi_quo_ih_4"], ax, c_values_b, 6)
axe(structure,deviation_enec_quo["wj_quo_ih_4"], ax, c_values_b, 7)
axe(structure,deviation_enec_quo["wk_quo_ih_4"], ax, c_values_b, 8)
axe(structure,deviation_enec_quo["wi_quo_jh_5"], ax, c_values_b, 9)
axe(structure,deviation_enec_quo["wj_quo_jh_5"], ax, c_values_b, 10)
axe(structure,deviation_enec_quo["wk_quo_jh_5"], ax, c_values_b, 11)
axe(structure,deviation_enec_quo["wi_quo_kh_6"], ax, c_values_b, 12)
axe(structure,deviation_enec_quo["wj_quo_kh_6"], ax, c_values_b, 13)
axe(structure,deviation_enec_quo["wk_quo_kh_6"], ax, c_values_b, 14)
axe(structure,deviation_enec_quo["wi_quo_F_7"], ax, c_values_b, 15)
axe(structure,deviation_enec_quo["wj_quo_F_7"], ax, c_values_b, 16)
axe(structure,deviation_enec_quo["wk_quo_F_7"], ax, c_values_b, 17)
axe(structure,ENEC_quo, ax, c_values_b, 18)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18])
ax.set_yticklabels(["$w_i$($\phi$-$f_{ij}$)", "$w_j$($\phi$-$f_{ij}$)","$w_i$($\phi$-$f_{ik}$)",
                    "$w_k$($\phi$-$f_{ik}$)","$w_j$($\phi$-$f_{jk}$)", "$w_k$($\phi$-$f_{jk}$)",
                    "$w_i$($\phi$-ih)", "$w_j$($\phi$-ih)", "$w_k$($\phi$-ih)", "$w_i$($\phi$-jh)", "$w_j$($\phi$-jh)",
                    "$w_k$($\phi$-jh)",
                    "$w_i$($\phi$-kh)", "$w_j$($\phi$-kh)", "$w_k$($\phi$-kh)", "$w_i$($\phi$-F)", "$w_j$($\phi$-F)",
                    "$w_k$($\phi$-F)", "ENEC_$\phi$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition bilatérale $\phi$")
plt.show()

# Déviation coalition - ih
# déviation manquante
dev_manquant_fij, dev_manquant_inverse_fij, dev_global_manquant_fij = deviation_uni("fij",
                            ['fik', 'fjk'], w_fij,[w_fik, w_fjk],  ['i','j','k'],
                                                                                    structure,  n_b, c_values_b)
dev_manquant_fik, dev_manquant_inverse_fik, dev_global_manquant_fik = deviation_uni("fik",
                            ['fij', 'fjk'], w_fik,[w_fij, w_fjk],  ['i','j','k'],
                                                                                    structure, n_b, c_values_b)

vect_dev_ih = [w_F, w_quo, w_jh, w_kh, w_fjk]
symbole_dev = ['F', 'quo', 'jh', 'kh', 'fjk']
vect_pays = [['j','k'], ['j','k'], ['j','k'],['j','k'],['j','k']]
vect_cre_ih = [[[dev_inverse_F['w_ih_F'][1], dev_inverse_F['w_kh_F'][1], dev_inverse_F['w_fik_F'][1]],
                [dev_inverse_F['w_ih_F'][2], dev_inverse_F['w_jh_F'][2], dev_inverse_F['w_fij_F'][2]]],
                [[dev_global_fij['w_fij_quo'][1]],[dev_global_fik['w_fik_quo'][2]]],
                [[dev_inverse_jh['w_fij_jh'][1], dev_inverse_jh['w_fjk_jh'][1], dev_inverse_jh['w_quo_jh'][1]],
                 [dev_inverse_jh['w_fij_jh'][2], dev_global_F['w_F_jh'][2]]],
                [[dev_inverse_kh['w_fik_kh'][1], dev_global_F['w_F_kh'][1]], [dev_inverse_kh['w_fik_kh'][2],
                                                        dev_inverse_kh['w_fjk_kh'][2], dev_inverse_kh['w_quo_kh'][2]]],
                [[dev_global_jh['w_jh_fjk'][1], dev_global_manquant_fij['w_fij_fjk'][1],
                  dev_inverse_fjk['w_quo_fjk'][1]],[dev_global_kh['w_kh_fjk'][2],
                                            dev_global_manquant_fik['w_fik_fjk'][2], dev_inverse_fjk['w_quo_fjk'][2]]]
                ]


ENEC_ih = []
deviation_enec_ih = deviation_coal('ih', symbole_dev, w_ih, vect_dev_ih, vect_cre_ih, vect_pays,
                                   structure, n_b, c_values_b)
for i in range(n_b):
    if EN_ih[i] == 0 :
        ENEC_ih.append(0)
    elif ((deviation_enec_ih['wj_ih_F_1'][i] + deviation_enec_ih['wk_ih_F_1'][i]) == 0
        or (deviation_enec_ih['wj_ih_quo_2'][i] + deviation_enec_ih['wk_ih_quo_2'][i]) == 1
        or (deviation_enec_ih['wj_ih_jh_3'][i] + deviation_enec_ih['wk_ih_jh_3'][i]) == 0
        or (deviation_enec_ih['wj_ih_kh_4'][i] + deviation_enec_ih['wk_ih_kh_4'][i]) == 0
        or (deviation_enec_ih['wj_ih_fjk_5'][i] + deviation_enec_ih['wk_ih_fjk_5'][i]) == 0):
        ENEC_ih.append(0)
    else: ENEC_ih.append(1)
print(ENEC_ih)
inverse_ENEC_ih = [1 - x for x in ENEC_ih]

fig, ax = plt.subplots()
axe(structure,deviation_enec_ih['wj_ih_F_1'], ax, c_values_b, 0)
axe(structure,deviation_enec_ih['wk_ih_F_1'], ax, c_values_b, 1)
axe(structure,deviation_enec_ih['wj_ih_quo_2'], ax, c_values_b, 2)
axe(structure,deviation_enec_ih['wk_ih_quo_2'], ax, c_values_b, 3)
axe(structure,deviation_enec_ih['wj_ih_jh_3'], ax, c_values_b, 4)
axe(structure,deviation_enec_ih['wk_ih_jh_3'], ax, c_values_b, 5)
axe(structure,deviation_enec_ih['wj_ih_kh_4'], ax, c_values_b, 6)
axe(structure,deviation_enec_ih['wk_ih_kh_4'], ax, c_values_b, 7)
axe(structure,deviation_enec_ih['wj_ih_fjk_5'], ax, c_values_b, 8)
axe(structure,deviation_enec_ih['wk_ih_fjk_5'], ax, c_values_b, 9)
axe(structure,ENEC_ih, ax, c_values_b, 10)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
ax.set_yticklabels(["$w_j$(ih-F)", "$w_k$(ih-F)", "$w_j$(ih-$\phi$)", "$w_k$(ih-$\phi$)", "$w_j$(ih-jh)",
                    "$w_k$(ih-jh)", "$w_j$(ih-kh)", "$w_k$(ih-kh)", "$w_j$(ih-$f_{jk}$)", "$w_k$(ih-$f_{jk}$)",
                    "ENEC_ih"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition bilatérale ih")
plt.show()

# Déviation coalition - jh
vect_dev_jh = [w_F, w_quo, w_ih, w_kh, w_fik]
symbole_dev = ['F', 'quo', 'ih', 'kh', 'fik']
vect_pays = [['i','k'], ['i','k'], ['i','k'],['i','k'],['i','k']]
vect_cre_jh = [[[dev_inverse_F['w_jh_F'][0], dev_inverse_F['w_kh_F'][0], dev_inverse_F['w_fjk_F'][0]],
                [dev_inverse_F['w_ih_F'][2], dev_inverse_F['w_jh_F'][2], dev_inverse_F['w_fij_F'][2]]],
                [[dev_global_fij['w_fij_quo'][0]],[dev_global_fjk['w_fjk_quo'][2]]],
                [[dev_inverse_ih['w_fij_ih'][0], dev_inverse_ih['w_fik_ih'][0], dev_inverse_ih['w_quo_ih'][0]],
                 [dev_inverse_ih['w_fij_ih'][2], dev_global_F['w_F_ih'][2]]],
                [[dev_inverse_kh['w_fjk_kh'][0], dev_global_F['w_F_kh'][0]], [dev_inverse_kh['w_fik_kh'][2],
                                                        dev_inverse_kh['w_fjk_kh'][2], dev_inverse_kh['w_quo_kh'][2]]],
                [[dev_global_ih['w_ih_fik'][0], dev_global_manquant_fij['w_fij_fik'][0],
                  dev_inverse_fik['w_quo_fik'][0]],[dev_global_kh['w_kh_fik'][2],
                                        dev_manquant_inverse_fik['w_fjk_fik'][2], dev_inverse_fik['w_quo_fik'][2]]]
                ]

ENEC_jh = []
deviation_enec_jh = deviation_coal('jh', symbole_dev, w_jh, vect_dev_jh, vect_cre_jh, vect_pays,
                                   structure, n_b, c_values_b)
for i in range(n_b):
    if EN_jh[i] == 0:
        ENEC_jh.append(0)
    elif ((deviation_enec_jh['wi_jh_F_1'][i] + deviation_enec_jh['wk_jh_F_1'][i]) == 0
        or (deviation_enec_jh['wi_jh_quo_2'][i] + deviation_enec_jh['wk_jh_quo_2'][i]) == 0
        or (deviation_enec_jh['wi_jh_ih_3'][i] + deviation_enec_jh['wk_jh_ih_3'][i]) == 0
        or (deviation_enec_jh['wi_jh_kh_4'][i] + deviation_enec_jh['wk_jh_kh_4'][i]) == 0
        or (deviation_enec_jh['wi_jh_fik_5'][i] + deviation_enec_jh['wk_jh_fik_5'][i]) == 0):
        ENEC_jh.append(0)
    else: ENEC_jh.append(1)
print(ENEC_jh)
inverse_ENEC_jh = [1 - x for x in ENEC_jh]

fig, ax = plt.subplots()
axe(structure,deviation_enec_jh['wi_jh_F_1'], ax, c_values_b, 0)
axe(structure,deviation_enec_jh['wk_jh_F_1'], ax, c_values_b, 1)
axe(structure,deviation_enec_jh['wi_jh_quo_2'], ax, c_values_b, 2)
axe(structure,deviation_enec_jh['wk_jh_quo_2'], ax, c_values_b, 3)
axe(structure,deviation_enec_jh['wi_jh_ih_3'], ax, c_values_b, 4)
axe(structure,deviation_enec_jh['wk_jh_ih_3'], ax, c_values_b, 5)
axe(structure,deviation_enec_jh['wi_jh_kh_4'], ax, c_values_b, 6)
axe(structure,deviation_enec_jh['wk_jh_kh_4'], ax, c_values_b, 7)
axe(structure,deviation_enec_jh['wi_jh_fik_5'], ax, c_values_b, 8)
axe(structure,deviation_enec_jh['wk_jh_fik_5'], ax, c_values_b, 9)
axe(structure,ENEC_jh, ax, c_values_b, 10)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
ax.set_yticklabels(["$w_i$(jh-F)", "$w_k$(jh-F)", "$w_i$(jh-$\phi$)", "$w_k$(jh-$\phi$)", "$w_i$(jh-ih)",
                    "$w_k$(jh-ih)", "$w_i$(jh-kh)", "$w_k$(jh-kh)", "$w_i$(jh-$f_{ik}$)", "$w_k$(jh-$f_{ik}$)",
                    "ENEC_jh"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition bilatérale jh")
plt.show()


# Déviation coalition - kh
vect_dev_kh = [w_F, w_quo, w_ih, w_jh, w_fij]
symbole_dev = ['F', 'quo', 'ih', 'jh', 'fij']
vect_pays = [['i','j'], ['i','j'], ['i','j'],['i','j'],['i','j']]
vect_cre_kh = [[[dev_inverse_F['w_kh_F'][0], dev_inverse_F['w_jh_F'][0], dev_inverse_F['w_fjk_F'][0]],
                [dev_inverse_F['w_ih_F'][1], dev_inverse_F['w_kh_F'][1], dev_inverse_F['w_fik_F'][1]]],
                [[dev_global_fik['w_fik_quo'][0]],[dev_global_fjk['w_fjk_quo'][1]]],
                [[dev_inverse_ih['w_fij_ih'][0], dev_inverse_ih['w_fik_ih'][0], dev_inverse_ih['w_quo_ih'][0]],
                 [dev_inverse_ih['w_fik_ih'][1], dev_global_F['w_F_ih'][1]]],
                [[dev_inverse_jh['w_fjk_jh'][0], dev_global_F['w_F_jh'][0]], [dev_inverse_jh['w_fij_jh'][1],
                                                        dev_inverse_jh['w_fjk_jh'][1], dev_inverse_jh['w_quo_jh'][1]]],
                [[dev_global_ih['w_ih_fij'][0], dev_global_manquant_fik['w_fik_fij'][0],
                  dev_inverse_fij['w_quo_fij'][0]],[dev_global_jh['w_jh_fij'][1],
                                        dev_manquant_inverse_fij['w_fjk_fij'][1], dev_inverse_fij['w_quo_fij'][1]]]
                ]

ENEC_kh = []
deviation_enec_kh = deviation_coal('kh', symbole_dev, w_kh, vect_dev_kh, vect_cre_kh, vect_pays,
                                   structure, n_b, c_values_b)
for i in range(n_b):
    if EN_kh[i] == 0 :
        ENEC_kh.append(0)
    elif ((deviation_enec_kh['wi_kh_F_1'][i] + deviation_enec_kh['wj_kh_F_1'][i]) == 0
        or (deviation_enec_kh['wi_kh_quo_2'][i] + deviation_enec_kh['wj_kh_quo_2'][i]) == 0
        or (deviation_enec_kh['wi_kh_ih_3'][i] + deviation_enec_kh['wj_kh_ih_3'][i]) == 0
        or (deviation_enec_kh['wi_kh_jh_4'][i] + deviation_enec_kh['wj_kh_jh_4'][i]) == 0
        or (deviation_enec_kh['wi_kh_fij_5'][i] + deviation_enec_kh['wj_kh_fij_5'][i]) == 0):
        ENEC_kh.append(0)
    else: ENEC_kh.append(1)
print(ENEC_kh)
inverse_ENEC_kh = [1 - x for x in ENEC_kh]

fig, ax = plt.subplots()
axe(structure,deviation_enec_kh['wi_kh_F_1'], ax, c_values_b, 0)
axe(structure,deviation_enec_kh['wj_kh_F_1'], ax, c_values_b, 1)
axe(structure,deviation_enec_kh['wi_kh_quo_2'], ax, c_values_b, 2)
axe(structure,deviation_enec_kh['wj_kh_quo_2'], ax, c_values_b, 3)
axe(structure,deviation_enec_kh['wi_kh_ih_3'], ax, c_values_b, 4)
axe(structure,deviation_enec_kh['wj_kh_ih_3'], ax, c_values_b, 5)
axe(structure,deviation_enec_kh['wi_kh_jh_4'], ax, c_values_b, 6)
axe(structure,deviation_enec_kh['wj_kh_jh_4'], ax, c_values_b, 7)
axe(structure,deviation_enec_kh['wi_kh_fij_5'], ax, c_values_b, 8)
axe(structure,deviation_enec_kh['wj_kh_fij_5'], ax, c_values_b, 9)
axe(structure,ENEC_kh, ax, c_values_b, 10)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
ax.set_yticklabels(["$w_i$(kh-F)", "$w_j$(kh-F)", "$w_i$(kh-$\phi$)", "$w_j$(kh-$\phi$)", "$w_i$(kh-ih)",
                "$w_j$(kh-ih)", "$w_i$(kh-jh)", "$w_j$(kh-jh)", "$w_i$(kh-$f_{ij}$)", "$w_j$(kh-$f_{ij}$)", "ENEC_kh"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition bilatérale kh")
plt.show()

# Déviation coalition - fij
vect_dev_fij = [w_F, w_ih, w_jh, w_kh, w_fik, w_fjk]
symbole_dev = ['F', 'ih', 'jh', 'kh', 'fik', 'fjk']
vect_pays = [['i', 'j', 'k'], ['i','k'], ['j','k'],['i', 'j', 'k'],['i','k'], ['j','k']]
vect_cre_fij = [[[inverse_ENEC_F],[inverse_ENEC_F], [inverse_ENEC_F]],
                [[dev_inverse_ih['w_fij_ih'][0], dev_inverse_ih['w_fik_ih'][0], dev_inverse_ih['w_quo_ih'][0]],
                 [dev_inverse_ih['w_fij_ih'][2]]],
                [[dev_inverse_jh['w_fij_jh'][1], dev_inverse_jh['w_fjk_jh'][1], dev_inverse_jh['w_quo_jh'][1]],
                 [dev_inverse_jh['w_fij_jh'][2]]],
                [[inverse_ENEC_kh],[inverse_ENEC_kh], [inverse_ENEC_kh]],
                [[dev_global_ih['w_ih_fik'][0], dev_global_manquant_fij['w_fij_fik'][0],
                  dev_inverse_fik['w_quo_fik'][0]], [dev_inverse_fik['w_quo_fik'][2]]],
                [[dev_global_jh['w_jh_fjk'][1], dev_global_manquant_fij['w_fij_fjk'][1],
                  dev_inverse_fjk['w_quo_fjk'][1]], [dev_inverse_fjk['w_quo_fjk'][2]]],
                ]

ENEC_fij = []
deviation_enec_fij = deviation_coal('fij', symbole_dev, w_fij, vect_dev_fij, vect_cre_fij, vect_pays,
                                    structure, n_b, c_values_b)
for i in range(n_b):
    if EN_fij[i] == 0 :
        ENEC_fij.append(0)
    elif ((deviation_enec_fij['wi_fij_F_1'][i] + deviation_enec_fij['wj_fij_F_1'][i] +
           deviation_enec_fij['wk_fij_F_1'][i]) == 0
        or (deviation_enec_fij['wi_fij_ih_2'][i] + deviation_enec_fij['wk_fij_ih_2'][i]) == 0
        or (deviation_enec_fij['wj_fij_jh_3'][i] + deviation_enec_fij['wk_fij_jh_3'][i]) == 0
        or (deviation_enec_fij['wi_fij_kh_4'][i] + deviation_enec_fij['wj_fij_kh_4'][i] +
            deviation_enec_fij['wk_fij_kh_4'][i]) == 0
        or (deviation_enec_fij['wi_fij_fik_5'][i] + deviation_enec_fij['wk_fij_fik_5'][i]) == 0
        or (deviation_enec_fij['wj_fij_fjk_6'][i] + deviation_enec_fij['wk_fij_fjk_6'][i]) == 0):
        ENEC_fij.append(0)
    else: ENEC_fij.append(1)
print(ENEC_fij)

fig, ax = plt.subplots()
axe(structure,deviation_enec_fij["wi_fij_F_1"], ax, c_values_b, 0)
axe(structure,deviation_enec_fij["wj_fij_F_1"], ax, c_values_b, 1)
axe(structure,deviation_enec_fij["wk_fij_F_1"], ax, c_values_b, 2)
axe(structure,deviation_enec_fij["wi_fij_ih_2"], ax, c_values_b, 3)
axe(structure,deviation_enec_fij["wk_fij_ih_2"], ax, c_values_b, 4)
axe(structure,deviation_enec_fij["wj_fij_jh_3"], ax, c_values_b, 5)
axe(structure,deviation_enec_fij["wk_fij_jh_3"], ax, c_values_b, 6)
axe(structure,deviation_enec_fij["wi_fij_kh_4"], ax, c_values_b, 7)
axe(structure,deviation_enec_fij["wj_fij_kh_4"], ax, c_values_b, 8)
axe(structure,deviation_enec_fij["wk_fij_kh_4"], ax, c_values_b, 9)
axe(structure,deviation_enec_fij["wi_fij_fik_5"], ax, c_values_b, 10)
axe(structure,deviation_enec_fij["wk_fij_fik_5"], ax, c_values_b, 11)
axe(structure,deviation_enec_fij["wj_fij_fjk_6"], ax, c_values_b, 12)
axe(structure,deviation_enec_fij["wk_fij_fjk_6"], ax, c_values_b, 13)
axe(structure,ENEC_fij, ax, c_values_b, 14)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14])
ax.set_yticklabels(["$w_i$($f_{ij}$-F)", "$w_j$($f_{ij}$-F)","$w_k$($f_{ij}$-F)", "$w_i$($f_{ij}$-ih)",
                    "$w_k$($f_{ij}$-ih)", "$w_j$($f_{ij}$-jh)",
                    "$w_k$($f_{ij}$-jh)", "$w_i$($f_{ij}$-kh)", "$w_j$($f_{ij}$-kh)", "$w_k$($f_{ij}$-kh)",
                    "$w_i$($f_{ij}$-$f_{ik}$)", "$w_k$($f_{ij}$-$f_{ik}$)",
                    "$w_j$($f_{ij}$-$f_{jk}$)", "$w_k$($f_{ij}$-$f_{jk}$)", "ENEC_$f_{ij}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition bilatérale $f_{ij}$")
plt.show()


# Déviation coalition - fik
vect_dev_fik = [w_F, w_ih, w_jh, w_kh, w_fij, w_fjk]
symbole_dev = ['F', 'ih', 'jh', 'kh', 'fij', 'fjk']
vect_pays = [['i', 'j', 'k'], ['i','j'], ['i', 'j', 'k'], ['j','k'], ['i','j'], ['j','k']]
vect_cre_fik = [[[inverse_ENEC_F],[inverse_ENEC_F], [inverse_ENEC_F]],
                [[dev_inverse_ih['w_fij_ih'][0], dev_inverse_ih['w_fik_ih'][0], dev_inverse_ih['w_quo_ih'][0]],
                 [dev_inverse_ih['w_fik_ih'][1]]],
                [[inverse_ENEC_jh],[inverse_ENEC_jh], [inverse_ENEC_jh]],
                [[dev_inverse_kh['w_fik_kh'][1]], [dev_inverse_kh['w_fik_kh'][2], dev_inverse_kh['w_fjk_kh'][2],
                                                   dev_inverse_kh['w_quo_kh'][2]]],
                [[dev_global_ih['w_ih_fij'][0], dev_global_manquant_fik['w_fik_fij'][0],
                  dev_inverse_fij['w_quo_fij'][0]], [dev_inverse_fij['w_quo_fij'][1]]],
                [[dev_inverse_fjk['w_quo_fjk'][1]], [dev_global_kh['w_kh_fjk'][2],
                                            dev_global_manquant_fik['w_fik_fjk'][2], dev_inverse_fjk['w_quo_fjk'][2]]],
                ]

ENEC_fik = []
deviation_enec_fik = deviation_coal('fik', symbole_dev, w_fik, vect_dev_fik, vect_cre_fik, vect_pays,
                                    structure, n_b, c_values_b)
for i in range(n_b):
    if EN_fik[i] == 0 :
        ENEC_fik.append(0)
    elif ((deviation_enec_fik['wi_fik_F_1'][i] + deviation_enec_fik['wj_fik_F_1'][i] +
           deviation_enec_fik['wk_fik_F_1'][i]) == 0
        or (deviation_enec_fik['wi_fik_ih_2'][i] + deviation_enec_fik['wj_fik_ih_2'][i]) == 0
        or (deviation_enec_fik['wi_fik_jh_3'][i] + deviation_enec_fik['wj_fik_jh_3'][i] +
            deviation_enec_fik['wk_fik_jh_3'][i]) == 0
        or (deviation_enec_fik['wj_fik_kh_4'][i] + deviation_enec_fik['wk_fik_kh_4'][i]) == 0
        or (deviation_enec_fik['wi_fik_fij_5'][i] + deviation_enec_fik['wj_fik_fij_5'][i]) == 0
        or (deviation_enec_fik['wj_fik_fjk_6'][i] + deviation_enec_fik['wk_fik_fjk_6'][i]) == 0):
        ENEC_fik.append(0)
    else: ENEC_fik.append(1)
print(ENEC_fik)

fig, ax = plt.subplots()
axe(structure,deviation_enec_fik["wi_fik_F_1"], ax, c_values_b, 0)
axe(structure,deviation_enec_fik["wj_fik_F_1"], ax, c_values_b, 1)
axe(structure,deviation_enec_fik["wk_fik_F_1"], ax, c_values_b, 2)
axe(structure,deviation_enec_fik["wi_fik_ih_2"], ax, c_values_b, 3)
axe(structure,deviation_enec_fik["wj_fik_ih_2"], ax, c_values_b, 4)
axe(structure,deviation_enec_fik["wi_fik_jh_3"], ax, c_values_b, 5)
axe(structure,deviation_enec_fik["wj_fik_jh_3"], ax, c_values_b, 6)
axe(structure,deviation_enec_fik["wk_fik_jh_3"], ax, c_values_b, 7)
axe(structure,deviation_enec_fik["wj_fik_kh_4"], ax, c_values_b, 8)
axe(structure,deviation_enec_fik["wk_fik_kh_4"], ax, c_values_b, 9)
axe(structure,deviation_enec_fik["wi_fik_fij_5"], ax, c_values_b, 10)
axe(structure,deviation_enec_fik["wj_fik_fij_5"], ax, c_values_b, 11)
axe(structure,deviation_enec_fik["wj_fik_fjk_6"], ax, c_values_b, 12)
axe(structure,deviation_enec_fik["wk_fik_fjk_6"], ax, c_values_b, 13)
axe(structure,ENEC_fik, ax, c_values_b, 14)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14])
ax.set_yticklabels(["$w_i$($f_{ik}$-F)", "$w_j$($f_{ik}$-F)","$w_k$($f_{ik}$-F)", "$w_i$($f_{ik}$-ih)",
                    "$w_j$($f_{ik}$-ih)", "$w_i$($f_{ik}$-jh)",
                    "$w_j$($f_{ik}$-jh)", "$w_k$($f_{ik}$-jh)", "$w_j$($f_{ik}$-kh)", "$w_k$($f_{ik}$-kh)",
                    "$w_i$($f_{ik}$-$f_{ij}$)", "$w_j$($f_{ik}$-$f_{ij}$)",
                    "$w_j$($f_{ik}$-$f_{jk}$)", "$w_k$($f_{ik}$-$f_{jk}$)", "ENEC_$f_{ik}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition bilatérale $f_{ik}$")
plt.show()


# Déviation coalition - fjk
vect_dev_fjk = [w_F, w_ih, w_jh, w_kh, w_fij, w_fik]
symbole_dev = ['F', 'ih', 'jh', 'kh', 'fij', 'fik']
vect_pays = [['i', 'j', 'k'], ['i', 'j', 'k'], ['i','j'], ['i','k'], ['i','j'], ['i','k']]
vect_cre_fjk = [[[inverse_ENEC_F],[inverse_ENEC_F], [inverse_ENEC_F]],
                [[inverse_ENEC_ih],[inverse_ENEC_ih], [inverse_ENEC_ih]],
                [[dev_inverse_jh['w_fjk_jh'][0]], [dev_inverse_jh['w_fij_jh'][1], dev_inverse_jh['w_fjk_jh'][1],
                                                   dev_inverse_jh['w_quo_jh'][1]]],
                [[dev_inverse_kh['w_fjk_kh'][0]], [dev_inverse_kh['w_fik_kh'][2], dev_inverse_kh['w_fjk_kh'][2],
                                                   dev_inverse_kh['w_quo_kh'][2]]],
                [[dev_inverse_fij['w_quo_fij'][0]], [dev_global_jh['w_jh_fij'][1],
                                            dev_manquant_inverse_fij['w_fjk_fij'][1], dev_inverse_fij['w_quo_fij'][1]]],
                [[dev_inverse_fik['w_quo_fik'][0]], [dev_global_kh['w_kh_fik'][2],
                                            dev_manquant_inverse_fik['w_fjk_fik'][2], dev_inverse_fik['w_quo_fik'][2]]],
                ]


ENEC_fjk = []
deviation_enec_fjk = deviation_coal('fjk', symbole_dev, w_fjk, vect_dev_fjk, vect_cre_fjk,
                                    vect_pays, structure, n_b, c_values_b)
for i in range(n_b):
    if EN_fjk[i] == 0:
        ENEC_fjk.append(0)
    elif ((deviation_enec_fjk['wi_fjk_F_1'][i] + deviation_enec_fjk['wj_fjk_F_1'][i] +
           deviation_enec_fjk['wk_fjk_F_1'][i]) == 0
        or (deviation_enec_fjk['wi_fjk_ih_2'][i] + deviation_enec_fjk['wj_fjk_ih_2'][i] +
            deviation_enec_fjk['wk_fjk_ih_2'][i]) == 0
        or (deviation_enec_fjk['wi_fjk_jh_3'][i] + deviation_enec_fjk['wj_fjk_jh_3'][i]) == 0
        or (deviation_enec_fjk['wi_fjk_kh_4'][i] + deviation_enec_fjk['wk_fjk_kh_4'][i]) == 0
        or (deviation_enec_fjk['wi_fjk_fij_5'][i] + deviation_enec_fjk['wj_fjk_fij_5'][i]) == 0
        or (deviation_enec_fjk['wi_fjk_fik_6'][i] + deviation_enec_fjk['wk_fjk_fik_6'][i]) == 0):
        ENEC_fjk.append(0)
    else: ENEC_fjk.append(1)
print(ENEC_fjk)

fig, ax = plt.subplots()
axe(structure,deviation_enec_fjk["wi_fjk_F_1"], ax, c_values_b, 0)
axe(structure,deviation_enec_fjk["wj_fjk_F_1"], ax, c_values_b, 1)
axe(structure,deviation_enec_fjk["wk_fjk_F_1"], ax, c_values_b, 2)
axe(structure,deviation_enec_fjk["wi_fjk_ih_2"], ax, c_values_b, 3)
axe(structure,deviation_enec_fjk["wj_fjk_ih_2"], ax, c_values_b, 4)
axe(structure,deviation_enec_fjk["wk_fjk_ih_2"], ax, c_values_b, 5)
axe(structure,deviation_enec_fjk["wi_fjk_jh_3"], ax, c_values_b, 6)
axe(structure,deviation_enec_fjk["wj_fjk_jh_3"], ax, c_values_b, 7)
axe(structure,deviation_enec_fjk["wi_fjk_kh_4"], ax, c_values_b, 8)
axe(structure,deviation_enec_fjk["wk_fjk_kh_4"], ax, c_values_b, 9)
axe(structure,deviation_enec_fjk["wi_fjk_fij_5"], ax, c_values_b, 10)
axe(structure,deviation_enec_fjk["wj_fjk_fij_5"], ax, c_values_b, 11)
axe(structure,deviation_enec_fjk["wi_fjk_fik_6"], ax, c_values_b, 12)
axe(structure,deviation_enec_fjk["wk_fjk_fik_6"], ax, c_values_b, 13)
axe(structure,ENEC_fjk, ax, c_values_b, 14)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14])
ax.set_yticklabels(["$w_i$($f_{jk}$-F)", "$w_j$($f_{jk}$-F)","$w_k$($f_{jk}$-F)", "$w_i$($f_{jk}$-ih)",
                    "$w_j$($f_{jk}$-ih)", "$w_k$($f_{jk}$-ih)",
                    "$w_i$($f_{jk}$-jh)", "$w_j$($f_{jk}$-jh)", "$w_i$($f_{jk}$-kh)", "$w_k$($f_{jk}$-kh)",
                    "$w_i$($f_{jk}$-$f_{ij}$)", "$w_j$($f_{jk}$-$f_{ij}$)",
                    "$w_i$($f_{jk}$-$f_{ik}$)", "$w_k$($f_{jk}$-$f_{ik}$)", "ENEC_$f_{jk}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition bilatérale $f_{jk}$")
plt.show()


# Graphique EN bilatérale
fig, ax = plt.subplots()
axe(structure,ENEC_quo, ax, c_values_b, 0)
axe(structure,ENEC_fij, ax, c_values_b, 1)
axe(structure,ENEC_fik, ax, c_values_b, 2)
axe(structure,ENEC_fjk, ax, c_values_b, 3)
axe(structure,ENEC_ih, ax, c_values_b, 4)
axe(structure,ENEC_jh, ax, c_values_b, 5)
axe(structure,ENEC_kh, ax, c_values_b, 6)
axe(structure,ENEC_F, ax, c_values_b, 7)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_yticklabels(["ENEC_$\phi$", "ENEC_$f_{ij}$", "ENEC_$f_{ik}$", "ENEC_$f_{jk}$", "ENEC_ih", "ENEC_jh",
                    "ENEC_kh", "ENEC_F"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("ENEC - Déviation de coalition bilatérale")

plt.show()






##################################### Déviation de coalition - multilatérale ##############################################
structure = "multilatérale"

# Déviation coalition - F
vect_dev_m_F = [w_quo, w_quo, w_quo, w_mij, w_mik, w_mjk]
symbole_dev = ['quo', 'quo', 'quo', 'mij', 'mik', 'mjk']
vect_pays = [['i','j'], ['i','k'], ['j','k'],['i', 'j'],['i','k'], ['j','k']]
vect_cre_m_F = [[[dev_global_mik['w_mik_quo'][0]],[dev_global_mjk['w_mjk_quo'][1]]],
              [[dev_global_mij['w_mij_quo'][0]],[dev_global_mjk['w_mjk_quo'][2]]],
              [[dev_global_mij['w_mij_quo'][1]],[dev_global_mik['w_mik_quo'][2]]],
              [[dev_inverse_mij['w_quo_mij'][0]], [dev_inverse_mij['w_quo_mij'][1]]],
              [[dev_inverse_mik['w_quo_mik'][0]], [dev_inverse_mik['w_quo_mik'][2]]],
              [[dev_inverse_mjk['w_quo_mjk'][1]], [dev_inverse_mjk['w_quo_mjk'][2]]],
              ]

ENEC_m_F = []
deviation_enec_m_F = deviation_coal("mF",symbole_dev, w_F, vect_dev_m_F, vect_cre_m_F, vect_pays,
                                    structure, n_m, c_values_m)
for i in range(n_m):
    if EN_m_F[i] == 0:
        ENEC_m_F.append(0)
    elif ((deviation_enec_m_F['wi_mF_quo_1'][i] + deviation_enec_m_F['wj_mF_quo_1'][i]) == 0
            or (deviation_enec_m_F['wi_mF_quo_2'][i] + deviation_enec_m_F['wk_mF_quo_2'][i]) == 0
            or (deviation_enec_m_F['wj_mF_quo_3'][i] + deviation_enec_m_F['wk_mF_quo_3'][i]) == 0
            or (deviation_enec_m_F['wi_mF_mij_4'][i] + deviation_enec_m_F['wj_mF_mij_4'][i]) == 0
            or (deviation_enec_m_F['wi_mF_mik_5'][i] + deviation_enec_m_F['wk_mF_mik_5'][i]) == 0
            or (deviation_enec_m_F['wj_mF_mjk_6'][i] + deviation_enec_m_F['wk_mF_mjk_6'][i]) == 0):
        ENEC_m_F.append(0)
    else: ENEC_m_F.append(1)
print(ENEC_m_F)
inverse_ENEC_m_F = [1 - x for x in ENEC_m_F]

fig, ax = plt.subplots()
axe(structure,deviation_enec_m_F["wi_mF_quo_1"], ax, c_values_m, 0)
axe(structure,deviation_enec_m_F["wj_mF_quo_1"], ax, c_values_m, 1)
axe(structure,deviation_enec_m_F["wi_mF_quo_2"], ax, c_values_m, 2)
axe(structure,deviation_enec_m_F["wk_mF_quo_2"], ax, c_values_m, 3)
axe(structure,deviation_enec_m_F["wj_mF_quo_3"], ax, c_values_m, 4)
axe(structure,deviation_enec_m_F["wk_mF_quo_3"], ax, c_values_m, 5)
axe(structure,deviation_enec_m_F["wi_mF_mij_4"], ax, c_values_m, 6)
axe(structure,deviation_enec_m_F["wj_mF_mij_4"], ax, c_values_m, 7)
axe(structure,deviation_enec_m_F["wi_mF_mik_5"], ax, c_values_m, 8)
axe(structure,deviation_enec_m_F["wk_mF_mik_5"], ax, c_values_m, 9)
axe(structure,deviation_enec_m_F["wj_mF_mjk_6"], ax, c_values_m, 10)
axe(structure,deviation_enec_m_F["wk_mF_mjk_6"], ax, c_values_m, 11)
axe(structure,ENEC_m_F, ax, c_values_m, 12)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6,7,8,9,10,11,12])
ax.set_yticklabels(["$w_i$(F-$\phi$)", "$w_j$(F-$\phi$)","$w_i$(F-$\phi$)", "$w_k$(F-$\phi$)", "$w_j$(F-$\phi$)",
                    "$w_k$(F-$\phi$)",
                    "$w_i$(F-$m_{ij}$)", "$w_j$(F-$m_{ij}$)", "$w_i$(F-$m_{ik}$)","$w_k$(F-$m_{ik}$)",
                    "$w_j$(F-$m_{jk}$)", "$w_k$(F-$m_{jk}$)", "ENEC_F"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition multilatérale F")
plt.show()


# Déviation coalition - statu quo
vect_dev_m_quo = [w_mij, w_mik, w_mjk, w_F]
symbole_dev = ['mij', 'mik', 'mjk', 'mF']
vect_pays = [['i','j'], ['i','k'], ['j','k'], ['i', 'j', 'k']]
vect_cre_m_quo = [[[dev_inverse_mij['w_quo_mij'][0]],[dev_inverse_mij['w_quo_mij'][1]]],
                [[dev_inverse_mik['w_quo_mik'][0]],[dev_inverse_mik['w_quo_mik'][2]]],
                [[dev_inverse_mjk['w_quo_mjk'][1]],[dev_inverse_mjk['w_quo_mjk'][2]]],
                [[inverse_ENEC_m_F], [inverse_ENEC_m_F],[inverse_ENEC_m_F]]
                ]

ENEC_m_quo = []
deviation_enec_m_quo = deviation_coal('quo', symbole_dev, w_quo, vect_dev_m_quo, vect_cre_m_quo,
                                      vect_pays, structure, n_m, c_values_m)
for i in range(n_m):
    if EN_m_quo[i] == 0:
        ENEC_m_quo.append(0)
    elif ((deviation_enec_m_quo['wi_quo_mij_1'][i] + deviation_enec_m_quo['wj_quo_mij_1'][i]) == 0
        or (deviation_enec_m_quo['wi_quo_mik_2'][i] + deviation_enec_m_quo['wk_quo_mik_2'][i]) == 0
        or (deviation_enec_m_quo['wj_quo_mjk_3'][i] + deviation_enec_m_quo['wk_quo_mjk_3'][i]) == 0
        or (deviation_enec_m_quo['wj_quo_mF_4'][i] + deviation_enec_m_quo['wk_quo_mF_4'][i] +
            deviation_enec_m_quo['wj_quo_mF_4'][i]) == 0):
        ENEC_m_quo.append(0)
    else: ENEC_m_quo.append(1)
print(ENEC_m_quo)

fig, ax = plt.subplots()
axe(structure,deviation_enec_m_quo["wi_quo_mij_1"], ax, c_values_m, 0)
axe(structure,deviation_enec_m_quo["wj_quo_mij_1"], ax, c_values_m, 1)
axe(structure,deviation_enec_m_quo["wi_quo_mik_2"], ax, c_values_m, 2)
axe(structure,deviation_enec_m_quo["wk_quo_mik_2"], ax, c_values_m, 3)
axe(structure,deviation_enec_m_quo["wj_quo_mjk_3"], ax, c_values_m, 4)
axe(structure,deviation_enec_m_quo["wk_quo_mjk_3"], ax, c_values_m, 5)
axe(structure,deviation_enec_m_quo["wi_quo_mF_4"], ax, c_values_m, 6)
axe(structure,deviation_enec_m_quo["wj_quo_mF_4"], ax, c_values_m, 7)
axe(structure,deviation_enec_m_quo["wk_quo_mF_4"], ax, c_values_m, 8)
axe(structure,ENEC_m_quo, ax, c_values_m, 9)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7,8,9])
ax.set_yticklabels(["$w_i$($\phi$-$m_{ij}$)", "$w_j$($\phi$-$m_{ij}$)","$w_i$($\phi$-$m_{ik}$)",
                    "$w_k$($\phi$-$m_{ik}$)","$w_j$($\phi$-$m_{jk}$)", "$w_k$($\phi$-$m_{jk}$)",
                    "$w_i$($\phi$-F)", "$w_j$($\phi$-F)", "$w_k$($\phi$-F)", "ENEC_$\phi$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition multilatérale $\phi$")
plt.show()

# Déviation coalition - mij
dev_manquant_mij, dev_manquant_inverse_mij, dev_global_manquant_mij = deviation_uni("mij",
                                                        ['mik', 'mjk'], w_mij,[w_mik, w_mjk],
                                                                ['i','j','k'], structure, n_m, c_values_m)
dev_manquant_mik, dev_manquant_inverse_mik, dev_global_manquant_mik = deviation_uni("mik",
                                                        ['mij', 'mjk'], w_mik,[w_mij, w_mjk],
                                                            ['i','j','k'], structure, n_m, c_values_m)

vect_dev_mij = [w_F, w_mik, w_mjk]
symbole_dev = ['mF', 'mik', 'mjk']
vect_pays = [['i', 'j', 'k'], ['i','k'], ['j','k']]
vect_cre_mij = [[[inverse_ENEC_m_F],[inverse_ENEC_m_F], [inverse_ENEC_m_F]],
                [[dev_global_manquant_mij['w_mij_mik'][0], dev_inverse_mik['w_quo_mik'][0]],
                 [dev_inverse_mik['w_quo_mik'][2]]],
                [[dev_global_manquant_mij['w_mij_mjk'][1], dev_inverse_mjk['w_quo_mjk'][1]],
                 [dev_inverse_mjk['w_quo_mjk'][2]]],
                ]
ENEC_mij = []
deviation_enec_mij = deviation_coal('mij', symbole_dev, w_mij, vect_dev_mij, vect_cre_mij, vect_pays,
                                    structure, n_m, c_values_m)
for i in range(n_m):
    if EN_mij[i] == 0:
        ENEC_mij.append(0)
    elif ((deviation_enec_mij['wi_mij_mF_1'][i] + deviation_enec_mij['wj_mij_mF_1'][i] +
           deviation_enec_mij['wk_mij_mF_1'][i]) == 0
        or (deviation_enec_mij['wi_mij_mik_2'][i] + deviation_enec_mij['wk_mij_mik_2'][i]) == 0
        or (deviation_enec_mij['wj_mij_mjk_3'][i] + deviation_enec_mij['wk_mij_mjk_3'][i]) == 0):
        ENEC_mij.append(0)
    else: ENEC_mij.append(1)
print(ENEC_mij)

fig, ax = plt.subplots()
axe(structure,deviation_enec_mij["wi_mij_mF_1"], ax, c_values_m, 0)
axe(structure,deviation_enec_mij["wj_mij_mF_1"], ax, c_values_m, 1)
axe(structure,deviation_enec_mij["wk_mij_mF_1"], ax, c_values_m, 2)
axe(structure,deviation_enec_mij["wi_mij_mik_2"], ax, c_values_m, 3)
axe(structure,deviation_enec_mij["wk_mij_mik_2"], ax, c_values_m, 4)
axe(structure,deviation_enec_mij["wj_mij_mjk_3"], ax, c_values_m, 5)
axe(structure,deviation_enec_mij["wk_mij_mjk_3"], ax, c_values_m, 6)
axe(structure,ENEC_mij, ax, c_values_m, 7)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_yticklabels(["$w_i$($m_{ij}$-F)", "$w_j$($m_{ij}$-F)","$w_k$($m_{ij}$-F)", "$w_i$($m_{ij}$-$m_{ik}$)",
                    "$w_k$($m_{ij}$-$m_{ik}$)",
                    "$w_j$($m_{ij}$-$m_{jk}$)", "$w_k$($m_{ij}$-$m_{jk}$)", "ENEC_$m_{ij}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition multilatérale $m_{ij}$")
plt.show()


# Déviation coalition - mik
vect_dev_mik = [w_F, w_mij, w_mjk]
symbole_dev = ['mF', 'mij', 'mjk']
vect_pays = [['i', 'j', 'k'], ['i','j'], ['j','k']]
vect_cre_mik = [[[inverse_ENEC_m_F],[inverse_ENEC_m_F], [inverse_ENEC_m_F]],
                [[dev_global_manquant_mik['w_mik_mij'][0], dev_inverse_mij['w_quo_mij'][0]],
                 [dev_inverse_mij['w_quo_mij'][1]]],
                [[dev_inverse_mjk['w_quo_mjk'][1]], [dev_global_manquant_mik['w_mik_mjk'][2],
                                                     dev_inverse_mjk['w_quo_mjk'][2]]],
                ]

ENEC_mik = []
deviation_enec_mik = deviation_coal('mik', symbole_dev, w_mik, vect_dev_mik, vect_cre_mik, vect_pays,
                                    structure, n_m, c_values_m)
for i in range(n_m):
    if EN_mik[i] == 0:
        ENEC_mik.append(0)
    elif ((deviation_enec_mik['wi_mik_mF_1'][i] + deviation_enec_mik['wj_mik_mF_1'][i] +
           deviation_enec_mik['wk_mik_mF_1'][i]) == 0
        or (deviation_enec_mik['wi_mik_mij_2'][i] + deviation_enec_mik['wj_mik_mij_2'][i]) == 0
        or (deviation_enec_mik['wj_mik_mjk_3'][i] + deviation_enec_mik['wk_mik_mjk_3'][i]) == 0):
        ENEC_mik.append(0)
    else: ENEC_mik.append(1)
print(ENEC_mik)

fig, ax = plt.subplots()
axe(structure,deviation_enec_mik["wi_mik_mF_1"], ax, c_values_m, 0)
axe(structure,deviation_enec_mik["wj_mik_mF_1"], ax, c_values_m, 1)
axe(structure,deviation_enec_mik["wk_mik_mF_1"], ax, c_values_m, 2)
axe(structure,deviation_enec_mik["wi_mik_mij_2"], ax, c_values_m, 3)
axe(structure,deviation_enec_mik["wj_mik_mij_2"], ax, c_values_m, 4)
axe(structure,deviation_enec_mik["wj_mik_mjk_3"], ax, c_values_m, 5)
axe(structure,deviation_enec_mik["wk_mik_mjk_3"], ax, c_values_m, 6)
axe(structure,ENEC_mik, ax, c_values_m, 7)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_yticklabels(["$w_i$($m_{ik}$-F)", "$w_j$($m_{ik}$-F)","$w_k$($m_{ik}$-F)", "$w_i$($m_{ik}$-$m_{ij}$)",
                    "$w_k$($m_{ik}$-$m_{ij}$)",
                    "$w_j$($m_{ik}$-$m_{jk}$)", "$w_k$($m_{ik}$-$m_{jk}$)", "ENEC_$m_{ik}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition multilatérale $m_{ik}$")
plt.show()

# Déviation coalition - mjk
vect_dev_mjk = [w_F, w_mij, w_mik]
symbole_dev = ['mF', 'mij', 'mik']
vect_pays = [['i', 'j', 'k'], ['i','j'], ['i','k']]
vect_cre_mjk = [[[inverse_ENEC_m_F],[inverse_ENEC_m_F], [inverse_ENEC_m_F]],
                [[dev_inverse_mij['w_quo_mij'][0]], [dev_manquant_inverse_mij['w_mjk_mij'][1],
                                                     dev_inverse_mij['w_quo_mij'][1]]],
                [[dev_inverse_mik['w_quo_mik'][0]], [dev_manquant_inverse_mik['w_mjk_mik'][2],
                                                     dev_inverse_mik['w_quo_mik'][2]]],
                ]

ENEC_mjk = []
deviation_enec_mjk = deviation_coal('mjk', symbole_dev, w_mjk, vect_dev_mjk, vect_cre_mjk, vect_pays,
                                    structure, n_m, c_values_m)
for i in range(n_m):
    if EN_mjk[i] == 0:
        ENEC_mjk.append(0)
    elif ((deviation_enec_mjk['wi_mjk_mF_1'][i] + deviation_enec_mjk['wj_mjk_mF_1'][i] +
           deviation_enec_mjk['wk_mjk_mF_1'][i]) == 0
        or (deviation_enec_mjk['wi_mjk_mij_2'][i] + deviation_enec_mjk['wj_mjk_mij_2'][i]) == 0
        or (deviation_enec_mjk['wi_mjk_mik_3'][i] + deviation_enec_mjk['wk_mjk_mik_3'][i]) == 0):
        ENEC_mjk.append(0)
    else: ENEC_mjk.append(1)
print(ENEC_mjk)

fig, ax = plt.subplots()
axe(structure,deviation_enec_mjk["wi_mjk_mF_1"], ax, c_values_m, 0)
axe(structure,deviation_enec_mjk["wj_mjk_mF_1"], ax, c_values_m, 1)
axe(structure,deviation_enec_mjk["wk_mjk_mF_1"], ax, c_values_m, 2)
axe(structure,deviation_enec_mjk["wi_mjk_mij_2"], ax, c_values_m, 3)
axe(structure,deviation_enec_mjk["wj_mjk_mij_2"], ax, c_values_m, 4)
axe(structure,deviation_enec_mjk["wi_mjk_mik_3"], ax, c_values_m, 5)
axe(structure,deviation_enec_mjk["wk_mjk_mik_3"], ax, c_values_m, 6)
axe(structure,ENEC_mjk, ax, c_values_m, 7)
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_yticklabels(["$w_i$($m_{jk}$-F)", "$w_j$($m_{jk}$-F)","$w_k$($m_{jk}$-F)", "$w_i$($m_{jk}$-$m_{ij}$)",
                    "$w_j$($m_{jk}$-$m_{ij}$)",
                    "$w_i$($m_{jk}$-$m_{ik}$)", "$w_k$($m_{jk}$-$m_{ik}$)", "ENEC_$m_{jk}$"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("Déviation de coalition multilatérale $m_{jk}$")
plt.show()

# Graphique EN multilatérale
fig, ax = plt.subplots()
axe(structure,ENEC_m_quo, ax, c_values_m, 0)
axe(structure,ENEC_mij, ax, c_values_m, 1)
axe(structure,ENEC_mik, ax, c_values_m, 2)
axe(structure,ENEC_mjk, ax, c_values_m, 3)
axe(structure,ENEC_m_F, ax, c_values_m, 4)
ax.set_yticks([0, 1, 2, 3, 4])
ax.set_yticklabels(["ENEC_$\phi$", "ENEC_$m_{ij}$", "ENEC_$m_{ik}$", "ENEC_$m_{jk}$", "ENEC_F"])
ax.set_xlabel("Valeurs de $\lambda$")
ax.set_title("ENEC - Déviation de coalition multilatérale")
plt.show()
