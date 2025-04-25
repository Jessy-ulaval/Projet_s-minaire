import sys

sys.path.append("/fonctions_projet.py")
from fonctions_projet import *

############################# Situation 1 ##################################

# Tarifs optimaux
s_t_o_b_3_vect, s_t_o_b_3_dict = tarifs_optimaux_bilateraux(1,1, 1,ci, cj, ci, cj, 0,0)

s_t_o_m_3_vect, s_t_o_m_3_dict = tarifs_optimaux_multilateraux(1,1, 1,ci, cj, ci, cj, 0,0)


# Prix d'équilibre
p_equ_vect, p_equ_dict = prix_equilibre(a, 1, 1, 1, s_t_o_b_3_dict["ti_quo"], s_t_o_b_3_dict["tj_quo"],
                                        s_t_o_b_3_dict["tk_quo"], s_t_o_b_3_dict["ti_fij"], s_t_o_b_3_dict["tj_fij"]
                                        , s_t_o_b_3_dict["ti_fik"], s_t_o_b_3_dict["tk_fik"], s_t_o_b_3_dict["tj_fjk"],
                                        s_t_o_b_3_dict["tk_fjk"], s_t_o_m_3_dict["ti_mij"]
                                        , s_t_o_m_3_dict["tj_mij"], s_t_o_m_3_dict["ti_mik"], s_t_o_m_3_dict["tk_mik"],
                                        s_t_o_m_3_dict["tj_mjk"], s_t_o_m_3_dict["tk_mjk"]
                                        , s_t_o_b_3_dict["tj_ih"], s_t_o_b_3_dict["tk_ih"], s_t_o_b_3_dict["ti_jh"],
                                        s_t_o_b_3_dict["tk_jh"], s_t_o_b_3_dict["ti_kh"], s_t_o_b_3_dict["tj_kh"]
                                        ,ci, cj, ci, cj, 0,0)
prix_equilibre_1(a, 1, 1, 1, s_t_o_b_3_dict["ti_quo"], s_t_o_b_3_dict["tj_quo"],
                                        s_t_o_b_3_dict["tk_quo"], s_t_o_b_3_dict["ti_fij"], s_t_o_b_3_dict["tj_fij"]
                                        , s_t_o_b_3_dict["ti_fik"], s_t_o_b_3_dict["tk_fik"], s_t_o_b_3_dict["tj_fjk"],
                                        s_t_o_b_3_dict["tk_fjk"], s_t_o_m_3_dict["ti_mij"]
                                        , s_t_o_m_3_dict["tj_mij"], s_t_o_m_3_dict["ti_mik"], s_t_o_m_3_dict["tk_mik"],
                                        s_t_o_m_3_dict["tj_mjk"], s_t_o_m_3_dict["tk_mjk"]
                                        , s_t_o_b_3_dict["tj_ih"], s_t_o_b_3_dict["tk_ih"], s_t_o_b_3_dict["ti_jh"],
                                        s_t_o_b_3_dict["tk_jh"], s_t_o_b_3_dict["ti_kh"], s_t_o_b_3_dict["tj_kh"]
                                        ,ci, cj, ci, cj, 0,0)

# Exportation
x_quo = export_i(a, 1, 1, 1, p_equ_dict["pi_i_quo"], p_equ_dict["pj_j_quo"], p_equ_dict["pk_k_quo"],
                 s_t_o_b_3_dict["ti_quo"], s_t_o_b_3_dict["ti_quo"], s_t_o_b_3_dict["tj_quo"], s_t_o_b_3_dict["tj_quo"],
                 s_t_o_b_3_dict["tk_quo"], s_t_o_b_3_dict["tk_quo"],ci, cj, ci, cj, 0,0)
x_fij = export_i(a, 1, 1, 1, p_equ_dict["pi_i_fij"], p_equ_dict["pj_j_fij"], p_equ_dict["pk_k_fij"], 0,
                 s_t_o_b_3_dict["ti_fij"], 0, s_t_o_b_3_dict["tj_fij"], s_t_o_b_3_dict["tk_quo"],
                 s_t_o_b_3_dict["tk_quo"],ci, cj, ci, cj, 0,0)
x_fik = export_i(a, 1, 1, 1, p_equ_dict["pi_i_fik"], p_equ_dict["pj_j_fik"], p_equ_dict["pk_k_fik"],
                 s_t_o_b_3_dict["ti_fik"], 0, s_t_o_b_3_dict["tj_quo"], s_t_o_b_3_dict["tj_quo"], 0,
                 s_t_o_b_3_dict["tk_fik"],ci, cj, ci, cj, 0,0)
x_fjk = export_i(a, 1, 1, 1, p_equ_dict["pi_i_fjk"], p_equ_dict["pj_j_fjk"], p_equ_dict["pk_k_fjk"],
                 s_t_o_b_3_dict["ti_quo"], s_t_o_b_3_dict["ti_quo"], s_t_o_b_3_dict["tj_fjk"], 0,
                 s_t_o_b_3_dict["tk_fjk"], 0,ci, cj, ci, cj, 0,0)
print("")
x_mij = export_i(a, 1, 1, 1, p_equ_dict["pi_i_mij"], p_equ_dict["pj_j_mij"], p_equ_dict["pk_k_mij"],
                 s_t_o_m_3_dict["ti_mij"], s_t_o_m_3_dict["ti_mij"], s_t_o_m_3_dict["tj_mij"], s_t_o_m_3_dict["tj_mij"],
                 s_t_o_m_3_dict["tk_quo"], s_t_o_m_3_dict["tk_quo"],ci, cj, ci, cj, 0,0)
x_mik = export_i(a, 1, 1, 1, p_equ_dict["pi_i_mik"], p_equ_dict["pj_j_mik"], p_equ_dict["pk_k_mik"],
                 s_t_o_m_3_dict["ti_mik"], s_t_o_m_3_dict["ti_mik"], s_t_o_m_3_dict["tj_quo"], s_t_o_m_3_dict["tj_quo"],
                 s_t_o_m_3_dict["tk_mik"], s_t_o_m_3_dict["tk_mik"],ci, cj, ci, cj, 0,0)
x_mjk = export_i(a, 1, 1, 1, p_equ_dict["pi_i_mjk"], p_equ_dict["pj_j_mjk"], p_equ_dict["pk_k_mjk"],
                 s_t_o_m_3_dict["ti_quo"], s_t_o_m_3_dict["ti_quo"], s_t_o_m_3_dict["tj_mjk"], s_t_o_m_3_dict["tj_mjk"],
                 s_t_o_m_3_dict["tk_mjk"], s_t_o_m_3_dict["tk_mjk"],ci, cj, ci, cj, 0,0)
print("")
x_ih = export_i(a, 1, 1, 1, p_equ_dict["pi_i_ih"], p_equ_dict["pj_j_ih"], p_equ_dict["pk_k_ih"], 0, 0,
                0,
                s_t_o_b_3_dict["tj_ih"], 0, s_t_o_b_3_dict["tk_ih"],ci, cj, ci, cj, 0,0)
x_jh = export_i(a, 1, 1, 1, p_equ_dict["pi_i_jh"], p_equ_dict["pj_j_jh"], p_equ_dict["pk_k_jh"], 0,
                s_t_o_b_3_dict["ti_jh"], 0, 0, s_t_o_b_3_dict["tk_jh"], 0,ci, cj, ci, cj, 0,0)
x_kh = export_i(a, 1, 1, 1, p_equ_dict["pi_i_kh"], p_equ_dict["pj_j_kh"], p_equ_dict["pk_k_kh"],
                s_t_o_b_3_dict["ti_kh"], 0, s_t_o_b_3_dict["tj_kh"], 0, 0, 0,ci, cj, ci, cj, 0,
                0)
x_F = export_i(a, 1, 1, 1, p_equ_dict["pi_i_F"], p_equ_dict["pj_j_F"], p_equ_dict["pk_k_F"], 0, 0,
               0, 0, 0, 0,ci, cj, ci, cj, 0,0)
X_b, X_m = export_final(x_quo, x_fij, x_fik, x_fjk, x_mij, x_mik, x_mjk, x_ih, x_jh, x_kh, x_F)
print(x_fjk)
print(x_mjk)

obtenir_bornes_restrictives_3(s_t_o_b_3_vect, s_t_o_b_3_dict, s_t_o_m_3_vect, s_t_o_m_3_dict, X_b, X_m, ci, cj)

# Bien-être optimaux - Structure bilaterale
print("##### Statu quo ####")
w_quo = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, s_t_o_b_3_dict["ti_quo"],
                   s_t_o_b_3_dict["tj_quo"],
                   s_t_o_b_3_dict["ti_quo"], s_t_o_b_3_dict["tk_quo"], s_t_o_b_3_dict["tj_quo"],
                   s_t_o_b_3_dict["tk_quo"],"quo")
print("")

print("#### fij ####")
w_fij = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, 0, 0, s_t_o_b_3_dict["ti_fij"],
                   s_t_o_b_3_dict["tk_quo"],
                   s_t_o_b_3_dict["tj_fij"], s_t_o_b_3_dict["tk_quo"], "fij")
print("")

print("##### fik ####")
w_fik = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, s_t_o_b_3_dict["ti_fik"],
                   s_t_o_b_3_dict["tj_quo"], 0, 0,
                   s_t_o_b_3_dict["tj_quo"], s_t_o_b_3_dict["tk_fik"], "fik")
print("")

print("##### fjk ####")
w_fjk = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, s_t_o_b_3_dict["ti_quo"],
                   s_t_o_b_3_dict["tj_fjk"],
                   s_t_o_b_3_dict["ti_quo"], s_t_o_b_3_dict["tk_fjk"], 0, 0, "fjk")
print("")

print("##### ih ####")
w_ih = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, 0, 0, 0, 0, s_t_o_b_3_dict["tj_ih"],
                  s_t_o_b_3_dict["tk_ih"], "ih")
print("")

print("##### jh ####")
w_jh = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, 0, 0, s_t_o_b_3_dict["ti_jh"],
                  s_t_o_b_3_dict["tk_jh"], 0, 0, "jh")
print("")

print("##### kh ####")
w_kh = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, s_t_o_b_3_dict["ti_kh"], s_t_o_b_3_dict["tj_kh"],
                  0, 0, 0, 0, "kh")
print("")

print("##### F ####")
w_F = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, 0, 0, 0, 0, 0, 0, "F")
print("")

# Bien-être optimaux - Structure multilateral
print("##### mij ####")
w_mij = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, s_t_o_m_3_dict["ti_mij"],
                   s_t_o_m_3_dict["tj_mij"],
                   s_t_o_m_3_dict["ti_mij"], s_t_o_m_3_dict["tk_quo"], s_t_o_m_3_dict["tj_mij"],
                   s_t_o_m_3_dict["tk_quo"], "mij")
print("")

print("##### mik ####")
w_mik = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0,s_t_o_m_3_dict["ti_mik"],s_t_o_m_3_dict["tj_quo"],
                   s_t_o_m_3_dict["ti_mik"], s_t_o_m_3_dict["tk_mik"], s_t_o_m_3_dict["tj_quo"],
                   s_t_o_m_3_dict["tk_mik"], "mik")
print("")

print("##### mjk ####")
w_mjk = w_optimaux(a, 1, 1, 1, ci, cj, ci, cj, 0,0, s_t_o_m_3_dict["ti_quo"],
                   s_t_o_m_3_dict["tj_mjk"],
                   s_t_o_m_3_dict["ti_quo"], s_t_o_m_3_dict["tk_mjk"], s_t_o_m_3_dict["tj_mjk"],
                   s_t_o_m_3_dict["tk_mjk"], "mjk")

######################################### Déviation unilatérale - Bilatérale ########################################
c_values_b_ci = np.arange(0, 1, 0.01)
c_values_b_cj = np.arange(0, 2/5, 0.01)
n_b_ci = len(c_values_b_ci)
n_b_cj = len(c_values_b_cj)
vect_pays = ['i','j','k']
titre_b = 'unilatérale bilatérale'
titre_m = 'unilatérale multilatérale'
structure = "bilatérale"


# Déviation unilatérale - fij
vect_dev_fij = [w_quo]
symbole_dev = ['quo']
EN_fij = np.zeros((n_b_ci, n_b_cj))
deviation_fij, dev_inverse_fij = deviation_uni_matrice("fij",symbole_dev, w_fij, vect_dev_fij,vect_pays,
                                                    'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            EN_fij[i, j] = 2
        elif deviation_fij['wi_fij_quo'][i,j] == 1 and deviation_fij['wj_fij_quo'][i,j] == 1:
            EN_fij[i,j] = 1
        else: EN_fij[i,j] = 0
print(EN_fij)
inverse_fij = 1- EN_fij


situation_initiale = 'fij'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_fij)
graph_matrice(structure,n_b_ci, n_b_cj, EN_fij, "EN_$f_{ij}$")





# Déviation unilatérale - fik
vect_dev_fik = [w_quo]
symbole_dev = ['quo']
EN_fik = np.zeros((n_b_ci, n_b_cj))
deviation_fik, dev_inverse_fik = deviation_uni_matrice("fik",symbole_dev, w_fik, vect_dev_fik,vect_pays,
                                                    'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            EN_fik[i, j] = 2
        elif deviation_fik['wi_fik_quo'][i,j] == 1 and deviation_fik['wk_fik_quo'][i,j] == 1:
            EN_fik[i,j] = 1
        else: EN_fik[i,j] = 0
print(EN_fik)
inverse_fik = 1 - EN_fik

situation_initiale = 'fik'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_fik)
graph_matrice(structure,n_b_ci, n_b_cj, EN_fik, "EN_$f_{ik}$")



# Déviation unilatérale - fjk
vect_dev_fjk = [w_quo]
symbole_dev = ['quo']
EN_fjk = np.zeros((n_b_ci, n_b_cj))
deviation_fjk, dev_inverse_fjk = deviation_uni_matrice("fjk",symbole_dev, w_fjk, vect_dev_fjk, vect_pays,
                                                'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            EN_fjk[i, j] = 2
        elif deviation_fjk['wj_fjk_quo'][i,j] == 1 and deviation_fjk['wk_fjk_quo'][i,j] == 1:
            EN_fjk[i,j] = 1
        else: EN_fjk[i,j] = 0
print(EN_fjk)
inverse_fjk = 1 - EN_fjk

situation_initiale = 'fjk'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_fjk)
graph_matrice(structure,n_b_ci, n_b_cj, EN_fjk, "EN_$f_{jk}$")



# Déviation unilatérale - ih
vect_dev_ih = [w_fij, w_fik, w_quo]
symbole_dev = ['fij', 'fik', 'quo']
EN_ih = np.zeros((n_b_ci, n_b_cj))
deviation_ih, dev_inverse_ih = deviation_uni_matrice("ih",symbole_dev, w_ih, vect_dev_ih, vect_pays,
                                            'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            EN_ih[i, j] = 2
        elif (deviation_ih['wi_ih_fij'][i,j] == 1 and deviation_ih['wi_ih_fik'][i,j] == 1 and
              deviation_ih['wi_ih_quo'][i,j] == 1
                and deviation_ih['wj_ih_fik'][i,j] == 1 and deviation_ih['wk_ih_fij'][i,j] == 1):
            EN_ih[i,j] = 1
        else: EN_ih[i,j] = 0
print(EN_ih)
inverse_ih = 1 - EN_ih

situation_initiale = 'ih'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_ih)
graph_matrice(structure,n_b_ci, n_b_cj, EN_ih, "EN_ih")


# Déviation unilatérale - jh
vect_dev_jh = [w_fjk, w_fij, w_quo]
symbole_dev = ['fjk', 'fij', 'quo']
EN_jh = np.zeros((n_b_ci, n_b_cj))
deviation_jh, dev_inverse_jh= deviation_uni_matrice("jh",symbole_dev, w_jh, vect_dev_jh, vect_pays,
                                        'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            EN_jh[i, j] = 2
        elif (deviation_jh['wj_jh_fjk'][i,j] == 1 and deviation_jh['wj_jh_fij'][i,j] == 1 and
              deviation_jh['wj_jh_quo'][i,j] == 1
                and deviation_jh['wk_jh_fij'][i,j] == 1 and deviation_jh['wi_jh_fjk'][i,j] == 1):
            EN_jh[i,j] = 1
        else: EN_jh[i,j] = 0
print(EN_jh)
inverse_jh = 1 - EN_jh

situation_initiale = 'jh'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_jh)
graph_matrice(structure,n_b_ci, n_b_cj, EN_jh, "EN_jh")



# Déviation unilatérale - kh
vect_dev_kh = [w_fjk, w_fik, w_quo]
symbole_dev = ['fjk', 'fik', 'quo']
EN_kh = np.zeros((n_b_ci, n_b_cj))
deviation_kh, dev_inverse_kh = deviation_uni_matrice("kh",symbole_dev, w_kh, vect_dev_kh, vect_pays,
                                    'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            EN_kh[i, j] = 2
        elif (deviation_kh['wk_kh_fjk'][i,j] == 1 and deviation_kh['wk_kh_fik'][i,j] == 1 and
              deviation_kh['wk_kh_quo'][i,j] == 1
                and deviation_kh['wj_kh_fik'][i,j] == 1 and deviation_kh['wi_kh_fjk'][i,j] == 1):
            EN_kh[i,j] = 1
        else: EN_kh[i,j] = 0
print(EN_kh)
inverse_kh = 1 - EN_kh

situation_initiale = 'kh'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_kh)
graph_matrice(structure,n_b_ci, n_b_cj, EN_kh, "EN_kh")


# Déviation unilatérale - F
vect_dev_F = [w_jh, w_kh, w_fjk, w_ih, w_fik, w_fij]
symbole_dev = ['jh', 'kh', 'fjk', 'ih', 'fik', 'fij']
EN_F = np.zeros((n_b_ci, n_b_cj))
deviation_F, dev_inverse_F = deviation_uni_matrice("F",symbole_dev, w_F, vect_dev_F, vect_pays,
                                    'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            EN_F[i, j] = 2
        elif (deviation_F['wi_F_jh'][i, j] == 1 and deviation_F['wi_F_kh'][i, j] == 1 and
              deviation_F['wi_F_fjk'][i, j] == 1 and
                deviation_F['wj_F_ih'][i, j] == 1 and deviation_F['wj_F_kh'][i, j] == 1 and
              deviation_F['wj_F_fik'][i, j] == 1 and
                deviation_F['wk_F_ih'][i, j] == 1 and deviation_F['wk_F_jh'][i, j] == 1 and
              deviation_F['wk_F_fij'][i, j] == 1):
            EN_F[i, j] = 1
        else: EN_F[i, j] = 0
print(EN_F)
inverse_F = 1 - EN_F

situation_initiale = 'F'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_F)
graph_matrice(structure,n_b_ci, n_b_cj, EN_F, "EN_F")





######################################### Déviation unilatérale - multilatérale #######################################
c_values_m_ci = np.arange(0, 1, 0.01)
c_values_m_cj = np.arange(0, 2/5, 0.01)
n_m_ci = len(c_values_m_ci)
n_m_cj = len(c_values_m_cj)
structure = "multilatérale"

# Déviation unilatérale - mij
vect_dev_mij = [w_quo]
symbole_dev = ['quo']
EN_mij = np.zeros((n_m_ci, n_m_cj))
deviation_mij, dev_inverse_mij = deviation_uni_matrice("mij",symbole_dev, w_mij, vect_dev_mij,vect_pays,
                                        'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
for i in range(n_m_ci):
    for j in range(n_m_cj):
        if deviation_mij['wi_mij_quo'][i, j] == 1 and deviation_mij['wj_mij_quo'][i, j] == 1:
            EN_mij[i, j] = 1
        else: EN_mij[i, j] = 0
print(EN_mij)
inverse_mij = 1 - EN_mij

situation_initiale = 'mij'
graph_matrice_e(structure,n_m_ci, n_m_cj, titre_m, situation_initiale,deviation_mij)
graph_matrice(structure,n_m_ci, n_m_cj, EN_mij, "EN_$m_{ij}$")

# Déviation unilatérale - mik
vect_dev_mik = [w_quo]
symbole_dev = ['quo']
EN_mik = np.zeros((n_m_ci, n_m_cj))
deviation_mik, dev_inverse_mik = deviation_uni_matrice("mik",symbole_dev, w_mik, vect_dev_mik,vect_pays,
                                    'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
for i in range(n_m_ci):
    for j in range(n_m_cj):
        if deviation_mik['wi_mik_quo'][i, j] == 1 and deviation_mik['wk_mik_quo'][i, j] == 1:
            EN_mik[i, j] = 1
        else: EN_mik[i, j] = 0
print(EN_mik)
inverse_mik = 1 - EN_mik

situation_initiale = 'mik'
graph_matrice_e(structure,n_m_ci, n_m_cj, titre_m, situation_initiale,deviation_mik)
graph_matrice(structure,n_m_ci, n_m_cj, EN_mik, "EN_$m_{ik}$")


# Déviation unilatérale - mjk
vect_dev_mjk = [w_quo]
symbole_dev = ['quo']
EN_mjk = np.zeros((n_m_ci, n_m_cj))
deviation_mjk, dev_inverse_mjk = deviation_uni_matrice("mjk",symbole_dev, w_mjk, vect_dev_mjk, vect_pays,
                                    'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
for i in range(n_m_ci):
    for j in range(n_m_cj):
        if deviation_mjk['wj_mjk_quo'][i, j] == 1 and deviation_mjk['wk_mjk_quo'][i, j] == 1:
            EN_mjk[i, j] = 1
        else: EN_mjk[i, j] = 0
print(EN_mjk)
inverse_mjk = 1 - EN_mjk

situation_initiale = 'mjk'
graph_matrice_e(structure,n_m_ci, n_m_cj, titre_m, situation_initiale,deviation_mjk)
graph_matrice(structure,n_m_ci, n_m_cj, EN_mjk, "EN_$m_{jk}$")



# Déviation unilatérale - F
vect_dev_F = [w_mjk, w_mik, w_mij]
symbole_dev = ['mjk', 'mik', 'mij']
EN_m_F = np.zeros((n_m_ci, n_m_cj))
deviation_m_F, dev_inverse_m_F = deviation_uni_matrice("F",symbole_dev, w_F, vect_dev_F, vect_pays,
                                        'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
for i in range(n_m_ci):
    for j in range(n_m_cj):
        if (deviation_m_F['wi_F_mjk'][i, j] == 1 and deviation_m_F['wj_F_mik'][i, j] == 1 and
            deviation_m_F['wk_F_mij'][i, j] == 1):
            EN_m_F[i, j] = 1
        else: EN_m_F[i, j] = 0
print(EN_m_F)
inverse_m_F = 1 - EN_m_F

situation_initiale = 'F'
graph_matrice_e(structure,n_m_ci, n_m_cj, titre_m, situation_initiale,deviation_m_F)
graph_matrice(structure,n_m_ci, n_m_cj, EN_m_F, "EN_F")



##################################### Déviation de coalition - Bilatérale ##############################################
titre_b = 'de coalition bilatérale'
titre_m = 'de coalition multilatérale'
structure = "bilatérale"
# Déviation coalition - F
vect_dev_F = [w_quo, w_quo, w_quo, w_fij, w_fik, w_fjk]
symbole_dev = ['quo', 'quo', 'quo', 'fij', 'fik', 'fjk']
vect_pays = [['i','j'], ['i','k'], ['j','k'],['i', 'j'],['i','k'], ['j','k']]
vect_cre_F = [[[deviation_fik['wi_fik_quo']],[deviation_fjk['wj_fjk_quo']]],
              [[deviation_fij['wi_fij_quo']],[deviation_fjk['wk_fjk_quo']]],
              [[deviation_fij['wj_fij_quo']],[deviation_fik['wk_fik_quo']]],
              [[deviation_ih['wi_ih_fij'],dev_inverse_fij['wi_quo_fij']], [deviation_jh['wj_jh_fij'],
                                                                           dev_inverse_fij['wj_quo_fij']]],
              [[deviation_ih['wi_ih_fik'],dev_inverse_fik['wi_quo_fik']], [deviation_kh['wk_kh_fik'],
                                                                           dev_inverse_fik['wk_quo_fik']]],
              [[deviation_jh['wj_jh_fjk'],dev_inverse_fjk['wj_quo_fjk']], [deviation_kh['wk_kh_fjk'],
                                                                           dev_inverse_fjk['wk_quo_fjk']]],
              ]

ENEC_F = np.zeros((n_b_ci, n_b_cj))
deviation_enec_F = deviation_coal_matrice("F",symbole_dev, w_F, vect_dev_F, vect_cre_F, vect_pays,
                                'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            ENEC_F[i, j] = 2
        elif EN_F[i,j] == 0:
            ENEC_F[i,j] = 0
        elif ((deviation_enec_F['wi_F_quo_1'][i, j] + deviation_enec_F['wj_F_quo_1'][i, j]) == 0
              or (deviation_enec_F['wi_F_quo_2'][i, j] + deviation_enec_F['wk_F_quo_2'][i, j]) == 0
              or (deviation_enec_F['wj_F_quo_3'][i, j] + deviation_enec_F['wk_F_quo_3'][i, j]) == 0
              or (deviation_enec_F['wi_F_fij_4'][i, j] + deviation_enec_F['wj_F_fij_4'][i, j]) == 0
              or (deviation_enec_F['wi_F_fik_5'][i, j] + deviation_enec_F['wk_F_fik_5'][i, j]) == 0
              or (deviation_enec_F['wj_F_fjk_6'][i, j] + deviation_enec_F['wk_F_fjk_6'][i, j]) == 0):
            ENEC_F[i, j] = 0
        else: ENEC_F[i, j] = 1
print(ENEC_F)
inverse_ENEC_F = 1 - ENEC_F

situation_initiale = 'F'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_enec_F)
graph_matrice(structure,n_b_ci, n_b_cj, ENEC_F, "ENEC_F")

deviation_enec_F_presentation = deviation_coal_matrice("F",symbole_dev, w_F, vect_dev_F, vect_cre_F,
                                                vect_pays, 'biilatérale', n_b_ci, n_b_cj, c_values_b_ci,
                                                       c_values_b_cj)
ENEC_F_final = np.zeros((n_b_ci, n_b_cj))
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            ENEC_F_final[i, j] = 0
        elif EN_F[i,j] == 0:
            ENEC_F_final[i,j] = 0
        elif ((deviation_enec_F_presentation['wi_F_quo_1'][i, j] +
               deviation_enec_F_presentation['wj_F_quo_1'][i, j]) == 0
              or (deviation_enec_F_presentation['wi_F_quo_2'][i, j] +
                  deviation_enec_F_presentation['wk_F_quo_2'][i, j]) == 0
              or (deviation_enec_F_presentation['wj_F_quo_3'][i, j] +
                  deviation_enec_F_presentation['wk_F_quo_3'][i, j]) == 0
              or (deviation_enec_F_presentation['wi_F_fij_4'][i, j] +
                  deviation_enec_F_presentation['wj_F_fij_4'][i, j]) == 0
              or (deviation_enec_F_presentation['wi_F_fik_5'][i, j] +
                  deviation_enec_F_presentation['wk_F_fik_5'][i, j]) == 0
              or (deviation_enec_F_presentation['wj_F_fjk_6'][i, j] +
                  deviation_enec_F_presentation['wk_F_fjk_6'][i, j]) == 0):
            ENEC_F_final[i, j] = 0
        else: ENEC_F_final[i, j] = 1
print(ENEC_F_final)


# Déviation coalition - statu quo
vect_dev_quo = [w_fij, w_fik, w_fjk, w_ih, w_jh, w_kh, w_F]
symbole_dev = ['fij', 'fik', 'fjk', 'ih', 'jh', 'kh', 'F']
vect_pays = [['i','j'], ['i','k'], ['j','k'],['i', 'j', 'k'],['i', 'j', 'k'], ['i', 'j', 'k'], ['i', 'j', 'k']]
vect_cre_quo = [[[dev_inverse_fij['wi_quo_fij']],[dev_inverse_fij['wj_quo_fij']]],
                [[dev_inverse_fik['wi_quo_fik']],[dev_inverse_fik['wk_quo_fik']]],
                [[dev_inverse_fjk['wj_quo_fjk']],[dev_inverse_fjk['wk_quo_fjk']]],
                [[inverse_ih], [inverse_ih, deviation_F['wj_F_ih']], [inverse_ih, deviation_F['wk_F_ih']]],
                [[inverse_jh, deviation_F['wi_F_jh']], [inverse_jh],  [inverse_jh, deviation_F['wk_F_jh']]],
                [[inverse_kh, deviation_F['wi_F_kh']], [inverse_kh, deviation_F['wj_F_kh']], [inverse_kh]],
                [[inverse_ENEC_F], [inverse_ENEC_F],[inverse_ENEC_F]]
                ]


ENEC_quo = np.zeros((n_b_ci, n_b_cj))
deviation_enec_quo = deviation_coal_matrice('quo', symbole_dev, w_quo, vect_dev_quo, vect_cre_quo,
                                            vect_pays, 'bilatérale', n_b_ci, n_b_cj, c_values_b_ci,
                                            c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            ENEC_quo[i, j] = 2
        elif ((deviation_enec_quo['wi_quo_fij_1'][i, j] + deviation_enec_quo['wj_quo_fij_1'][i, j]) == 0
              or (deviation_enec_quo['wi_quo_fik_2'][i, j] + deviation_enec_quo['wk_quo_fik_2'][i, j]) == 0
              or (deviation_enec_quo['wj_quo_fjk_3'][i, j] + deviation_enec_quo['wk_quo_fjk_3'][i, j]) == 0
              or (deviation_enec_quo['wi_quo_ih_4'][i, j] + deviation_enec_quo['wj_quo_ih_4'][i, j] +
                  deviation_enec_quo['wk_quo_ih_4'][i, j]) == 0
              or (deviation_enec_quo['wi_quo_jh_5'][i, j] + deviation_enec_quo['wj_quo_jh_5'][i, j] +
                  deviation_enec_quo['wk_quo_jh_5'][i, j]) == 0
              or (deviation_enec_quo['wi_quo_kh_6'][i, j] + deviation_enec_quo['wj_quo_kh_6'][i, j] +
                  deviation_enec_quo['wk_quo_kh_6'][i, j]) == 0
              or (deviation_enec_quo['wi_quo_F_7'][i, j] + deviation_enec_quo['wj_quo_F_7'][i, j] +
                  deviation_enec_quo['wk_quo_F_7'][i, j]) == 0):
            ENEC_quo[i, j] = 0
        else: ENEC_quo[i, j] = 1
print(ENEC_quo)

situation_initiale = 'statu quo'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_enec_quo)
graph_matrice(structure,n_b_ci, n_b_cj, ENEC_quo, "ENEC_$\phi$")


# Déviation coalition - ih
# déviation manquante
dev_manquant_fij, dev_manquant_inverse_fij = deviation_uni_matrice("fij",['fik', 'fjk'],
                                                w_fij,[w_fik, w_fjk],  ['i','j','k'],
                                                'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
dev_manquant_fik, dev_manquant_inverse_fik = deviation_uni_matrice("fik",['fij', 'fjk'],
                                                w_fik,[w_fij, w_fjk],  ['i','j','k'],
                                                    'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)

vect_dev_ih = [w_F, w_quo, w_jh, w_kh, w_fjk]
symbole_dev = ['F', 'quo', 'jh', 'kh', 'fjk']
vect_pays = [['j','k'], ['j','k'], ['j','k'],['j','k'],['j','k']]
vect_cre_ih = [[[dev_inverse_F['wj_ih_F'], dev_inverse_F['wj_kh_F'], dev_inverse_F['wj_fik_F']],
                [dev_inverse_F['wk_ih_F'], dev_inverse_F['wk_jh_F'], dev_inverse_F['wk_fij_F']]],
                [[deviation_fij['wj_fij_quo']],[deviation_fik['wk_fik_quo']]],
                [[dev_inverse_jh['wj_fij_jh'], dev_inverse_jh['wj_fjk_jh'], dev_inverse_jh['wj_quo_jh']],
                 [dev_inverse_jh['wk_fij_jh'], deviation_F['wk_F_jh']]],
                [[dev_inverse_kh['wj_fik_kh'], deviation_F['wj_F_kh']], [dev_inverse_kh['wk_fik_kh'],
                                                                         dev_inverse_kh['wk_fjk_kh'], dev_inverse_kh['wk_quo_kh']]],
                [[deviation_jh['wj_jh_fjk'], dev_manquant_fij['wj_fij_fjk'], dev_inverse_fjk['wj_quo_fjk']],
                 [deviation_kh['wk_kh_fjk'], dev_manquant_fik['wk_fik_fjk'], dev_inverse_fjk['wk_quo_fjk']]]
                ]

ENEC_ih = np.zeros((n_b_ci, n_b_cj))
deviation_enec_ih = deviation_coal_matrice('ih', symbole_dev, w_ih, vect_dev_ih, vect_cre_ih, vect_pays,
                            'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            ENEC_ih[i, j] = 2
        elif EN_ih[i, j] == 0:
            ENEC_ih[i, j] = 0
        elif ((deviation_enec_ih['wj_ih_F_1'][i, j] + deviation_enec_ih['wk_ih_F_1'][i, j]) == 0
              or (deviation_enec_ih['wj_ih_quo_2'][i, j] + deviation_enec_ih['wk_ih_quo_2'][i, j]) == 0
              or (deviation_enec_ih['wj_ih_jh_3'][i, j] + deviation_enec_ih['wk_ih_jh_3'][i, j]) == 0
              or (deviation_enec_ih['wj_ih_kh_4'][i, j] + deviation_enec_ih['wk_ih_kh_4'][i, j]) == 0
              or (deviation_enec_ih['wj_ih_fjk_5'][i, j] + deviation_enec_ih['wk_ih_fjk_5'][i, j]) == 0):
            ENEC_ih[i, j] = 0
        else: ENEC_ih[i, j] = 1
print(ENEC_ih)
inverse_ENEC_ih = 1 - ENEC_ih

situation_initiale = 'ih'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_enec_ih)
graph_matrice(structure,n_b_ci, n_b_cj, ENEC_ih, "ENEC_ih")




# Déviation coalition - jh
vect_dev_jh = [w_F, w_quo, w_ih, w_kh, w_fik]
symbole_dev = ['F', 'quo', 'ih', 'kh', 'fik']
vect_pays = [['i','k'], ['i','k'], ['i','k'],['i','k'],['i','k']]
vect_cre_jh = [[[dev_inverse_F['wi_jh_F'], dev_inverse_F['wi_kh_F'], dev_inverse_F['wi_fjk_F']],
                [dev_inverse_F['wk_ih_F'], dev_inverse_F['wk_jh_F'], dev_inverse_F['wk_fij_F']]],
                [[deviation_fij['wi_fij_quo']],[deviation_fjk['wk_fjk_quo']]],
                [[dev_inverse_ih['wi_fij_ih'], dev_inverse_ih['wi_fik_ih'], dev_inverse_ih['wi_quo_ih']],
                 [dev_inverse_ih['wk_fij_ih'], deviation_F['wk_F_ih']]],
                [[dev_inverse_kh['wi_fjk_kh'], deviation_F['wi_F_kh']], [dev_inverse_kh['wk_fik_kh'],
                                                                         dev_inverse_kh['wk_fjk_kh'], dev_inverse_kh['wk_quo_kh']]],
                [[deviation_ih['wi_ih_fik'], dev_manquant_fij['wi_fij_fik'], dev_inverse_fik['wi_quo_fik']],
                 [deviation_kh['wk_kh_fik'], dev_manquant_inverse_fik['wk_fjk_fik'], dev_inverse_fik['wk_quo_fik']]]
                ]

ENEC_jh = np.zeros((n_b_ci, n_b_cj))
deviation_enec_jh = deviation_coal_matrice('jh', symbole_dev, w_jh, vect_dev_jh, vect_cre_jh, vect_pays,
                                           'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            ENEC_jh[i, j] = 2
        elif EN_jh[i, j] == 0:
            ENEC_jh[i, j] = 0
        elif ((deviation_enec_jh['wi_jh_F_1'][i, j] + deviation_enec_jh['wk_jh_F_1'][i, j]) == 0
              or (deviation_enec_jh['wi_jh_quo_2'][i, j] + deviation_enec_jh['wk_jh_quo_2'][i, j]) == 0
              or (deviation_enec_jh['wi_jh_ih_3'][i, j] + deviation_enec_jh['wk_jh_ih_3'][i, j]) == 0
              or (deviation_enec_jh['wi_jh_kh_4'][i, j] + deviation_enec_jh['wk_jh_kh_4'][i, j]) == 0
              or (deviation_enec_jh['wi_jh_fik_5'][i, j] + deviation_enec_jh['wk_jh_fik_5'][i, j]) == 0):
            ENEC_jh[i, j] = 0
        else: ENEC_jh[i, j] = 1
print(ENEC_jh)
inverse_ENEC_jh = 1 - ENEC_jh

situation_initiale = 'jh'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_enec_jh)
graph_matrice(structure,n_b_ci, n_b_cj, ENEC_jh, "ENEC_jh")



# Déviation coalition - kh
vect_dev_kh = [w_F, w_quo, w_ih, w_jh, w_fij]
symbole_dev = ['F', 'quo', 'ih', 'jh', 'fij']
vect_pays = [['i','j'], ['i','j'], ['i','j'],['i','j'],['i','j']]
vect_cre_kh = [[[dev_inverse_F['wi_kh_F'], dev_inverse_F['wi_jh_F'], dev_inverse_F['wi_fjk_F']],
                [dev_inverse_F['wj_ih_F'], dev_inverse_F['wj_kh_F'], dev_inverse_F['wj_fik_F']]],
                [[deviation_fik['wi_fik_quo']],[deviation_fjk['wj_fjk_quo']]],
                [[dev_inverse_ih['wi_fij_ih'], dev_inverse_ih['wi_fik_ih'], dev_inverse_ih['wi_quo_ih']],
                 [dev_inverse_ih['wj_fik_ih'], deviation_F['wj_F_ih']]],
                [[dev_inverse_jh['wi_fjk_jh'], deviation_F['wi_F_jh']], [dev_inverse_jh['wj_fij_jh'],
                                                        dev_inverse_jh['wj_fjk_jh'], dev_inverse_jh['wj_quo_jh']]],
                [[deviation_ih['wi_ih_fij'], dev_manquant_fik['wi_fik_fij'], dev_inverse_fij['wi_quo_fij']],
                 [deviation_jh['wj_jh_fij'], dev_manquant_inverse_fij['wj_fjk_fij'], dev_inverse_fij['wj_quo_fij']]]
                ]

ENEC_kh = np.zeros((n_b_ci, n_b_cj))
deviation_enec_kh = deviation_coal_matrice('kh', symbole_dev, w_kh, vect_dev_kh, vect_cre_kh, vect_pays,
                            'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            ENEC_kh[i, j] = 2
        elif EN_kh[i, j] == 0:
            ENEC_kh[i, j] = 0
        elif ((deviation_enec_kh['wi_kh_F_1'][i, j] + deviation_enec_kh['wj_kh_F_1'][i, j]) == 0
              or (deviation_enec_kh['wi_kh_quo_2'][i, j] + deviation_enec_kh['wj_kh_quo_2'][i, j]) == 0
              or (deviation_enec_kh['wi_kh_ih_3'][i, j] + deviation_enec_kh['wj_kh_ih_3'][i, j]) == 0
              or (deviation_enec_kh['wi_kh_jh_4'][i, j] + deviation_enec_kh['wj_kh_jh_4'][i, j]) == 0
              or (deviation_enec_kh['wi_kh_fij_5'][i, j] + deviation_enec_kh['wj_kh_fij_5'][i, j]) == 0):
            ENEC_kh[i, j] = 0
        else: ENEC_kh[i, j] = 1
print(ENEC_kh)
inverse_ENEC_kh = 1 - ENEC_kh

situation_initiale = 'kh'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_enec_kh)
graph_matrice(structure,n_b_ci, n_b_cj, ENEC_kh, "ENEC_kh")



# Déviation coalition - fij
vect_dev_fij = [w_F, w_ih, w_jh, w_kh, w_fik, w_fjk]
symbole_dev = ['F', 'ih', 'jh', 'kh', 'fik', 'fjk']
vect_pays = [['i', 'j', 'k'], ['i','k'], ['j','k'],['i', 'j', 'k'],['i','k'], ['j','k']]
vect_cre_fij = [[[inverse_ENEC_F],[inverse_ENEC_F], [inverse_ENEC_F]],
                [[dev_inverse_ih['wi_fij_ih'], dev_inverse_ih['wi_fik_ih'], dev_inverse_ih['wi_quo_ih']],
                 [dev_inverse_ih['wk_fij_ih']]],
                [[dev_inverse_jh['wj_fij_jh'], dev_inverse_jh['wj_fjk_jh'], dev_inverse_jh['wj_quo_jh']],
                 [dev_inverse_jh['wk_fij_jh']]],
                [[inverse_ENEC_kh],[inverse_ENEC_kh], [inverse_ENEC_kh]],
                [[deviation_ih['wi_ih_fik'], dev_manquant_fij['wi_fij_fik'], dev_inverse_fik['wi_quo_fik']],
                 [dev_inverse_fik['wk_quo_fik']]],
                [[deviation_jh['wj_jh_fjk'], dev_manquant_fij['wj_fij_fjk'], dev_inverse_fjk['wj_quo_fjk']],
                 [dev_inverse_fjk['wk_quo_fjk']]],
                ]

ENEC_fij = np.zeros((n_b_ci, n_b_cj))
deviation_enec_fij = deviation_coal_matrice('fij', symbole_dev, w_fij, vect_dev_fij, vect_cre_fij,
                            vect_pays, 'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            ENEC_fij[i, j] = 2
        elif EN_fij[i, j] == 0:
            ENEC_fij[i, j] = 0
        elif ((deviation_enec_fij['wi_fij_F_1'][i, j] + deviation_enec_fij['wj_fij_F_1'][i, j] +
               deviation_enec_fij['wk_fij_F_1'][i, j]) == 0
              or (deviation_enec_fij['wi_fij_ih_2'][i, j] + deviation_enec_fij['wk_fij_ih_2'][i, j]) == 0
              or (deviation_enec_fij['wj_fij_jh_3'][i, j] + deviation_enec_fij['wk_fij_jh_3'][i, j]) == 0
              or (deviation_enec_fij['wi_fij_kh_4'][i, j] + deviation_enec_fij['wj_fij_kh_4'][i, j] +
                  deviation_enec_fij['wk_fij_kh_4'][i, j]) == 0
              or (deviation_enec_fij['wi_fij_fik_5'][i, j] + deviation_enec_fij['wk_fij_fik_5'][i, j]) == 0
              or (deviation_enec_fij['wj_fij_fjk_6'][i, j] + deviation_enec_fij['wk_fij_fjk_6'][i, j]) == 0):
            ENEC_fij[i, j] = 0
        else: ENEC_fij[i, j] = 1
print(ENEC_fij)

situation_initiale = 'fij'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_enec_fij)
graph_matrice(structure,n_b_ci, n_b_cj, ENEC_fij, "ENEC_$f_{ij}$")



# Déviation coalition - fik
vect_dev_fik = [w_F, w_ih, w_jh, w_kh, w_fij, w_fjk]
symbole_dev = ['F', 'ih', 'jh', 'kh', 'fij', 'fjk']
vect_pays = [['i', 'j', 'k'], ['i','j'], ['i', 'j', 'k'], ['j','k'], ['i','j'], ['j','k']]
vect_cre_fik = [[[inverse_ENEC_F],[inverse_ENEC_F], [inverse_ENEC_F]],
                [[dev_inverse_ih['wi_fij_ih'], dev_inverse_ih['wi_fik_ih'], dev_inverse_ih['wi_quo_ih']],
                 [dev_inverse_ih['wj_fik_ih']]],
                [[inverse_ENEC_jh],[inverse_ENEC_jh], [inverse_ENEC_jh]],
                [[dev_inverse_kh['wj_fik_kh']], [dev_inverse_kh['wk_fik_kh'], dev_inverse_kh['wk_fjk_kh'],
                                                 dev_inverse_kh['wk_quo_kh']]],
                [[deviation_ih['wi_ih_fij'], dev_manquant_fik['wi_fik_fij'], dev_inverse_fij['wi_quo_fij']],
                 [dev_inverse_fij['wj_quo_fij']]],
                [[dev_inverse_fjk['wj_quo_fjk']], [deviation_kh['wk_kh_fjk'], dev_manquant_fik['wk_fik_fjk'],
                                                   dev_inverse_fjk['wk_quo_fjk']]],
                ]

ENEC_fik = np.zeros((n_b_ci, n_b_cj))
deviation_enec_fik = deviation_coal_matrice('fik', symbole_dev, w_fik, vect_dev_fik, vect_cre_fik,
                            vect_pays, 'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            ENEC_fik[i, j] = 2
        elif EN_fik[i, j] == 0:
            ENEC_fik[i, j] = 0
        elif ((deviation_enec_fik['wi_fik_F_1'][i, j] + deviation_enec_fik['wj_fik_F_1'][i, j] +
               deviation_enec_fik['wk_fik_F_1'][i, j]) == 0
              or (deviation_enec_fik['wi_fik_ih_2'][i, j] + deviation_enec_fik['wj_fik_ih_2'][i, j]) == 0
              or (deviation_enec_fik['wi_fik_jh_3'][i, j] + deviation_enec_fik['wj_fik_jh_3'][i, j] +
                  deviation_enec_fik['wk_fik_jh_3'][i, j]) == 0
              or (deviation_enec_fik['wj_fik_kh_4'][i, j] + deviation_enec_fik['wk_fik_kh_4'][i, j]) == 0
              or (deviation_enec_fik['wi_fik_fij_5'][i, j] + deviation_enec_fik['wj_fik_fij_5'][i, j]) == 0
              or (deviation_enec_fik['wj_fik_fjk_6'][i, j] + deviation_enec_fik['wk_fik_fjk_6'][i, j]) == 0):
            ENEC_fik[i, j] = 0
        else: ENEC_fik[i, j] = 1
print(ENEC_fik)

situation_initiale = 'fik'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_enec_fik)
graph_matrice(structure,n_b_ci, n_b_cj, ENEC_fik, "ENEC_$f_{ik}$")


# Déviation coalition - fjk
vect_dev_fjk = [w_F, w_ih, w_jh, w_kh, w_fij, w_fik]
symbole_dev = ['F', 'ih', 'jh', 'kh', 'fij', 'fik']
vect_pays = [['i', 'j', 'k'], ['i', 'j', 'k'], ['i','j'], ['i','k'], ['i','j'], ['i','k']]
vect_cre_fjk = [[[inverse_ENEC_F],[inverse_ENEC_F], [inverse_ENEC_F]],
                [[inverse_ENEC_ih],[inverse_ENEC_ih], [inverse_ENEC_ih]],
                [[dev_inverse_jh['wi_fjk_jh']], [dev_inverse_jh['wj_fij_jh'], dev_inverse_jh['wj_fjk_jh'],
                                                 dev_inverse_jh['wj_quo_jh']]],
                [[dev_inverse_kh['wi_fjk_kh']], [dev_inverse_kh['wk_fik_kh'], dev_inverse_kh['wk_fjk_kh'],
                                                 dev_inverse_kh['wk_quo_kh']]],
                [[dev_inverse_fij['wi_quo_fij']], [deviation_jh['wj_jh_fij'], dev_manquant_inverse_fij['wj_fjk_fij'],
                                                   dev_inverse_fij['wj_quo_fij']]],
                [[dev_inverse_fik['wi_quo_fik']], [deviation_kh['wk_kh_fik'], dev_manquant_inverse_fik['wk_fjk_fik'],
                                                   dev_inverse_fik['wk_quo_fik']]],
                ]


ENEC_fjk = np.zeros((n_b_ci, n_b_cj))
deviation_enec_fjk = deviation_coal_matrice('fjk', symbole_dev, w_fjk, vect_dev_fjk, vect_cre_fjk,
                                vect_pays, 'bilatérale', n_b_ci, n_b_cj, c_values_b_ci, c_values_b_cj)
for i in range(n_b_ci):
    for j in range(n_b_cj):
        if j > 20:
            ENEC_fjk[i, j] = 2
        elif EN_fjk[i, j] == 0:
            ENEC_fjk[i, j] = 0
        elif ((deviation_enec_fjk['wi_fjk_F_1'][i, j] + deviation_enec_fjk['wj_fjk_F_1'][i, j] +
               deviation_enec_fjk['wk_fjk_F_1'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_ih_2'][i, j] + deviation_enec_fjk['wj_fjk_ih_2'][i, j] +
                  deviation_enec_fjk['wk_fjk_ih_2'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_jh_3'][i, j] + deviation_enec_fjk['wj_fjk_jh_3'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_kh_4'][i, j] + deviation_enec_fjk['wk_fjk_kh_4'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_fij_5'][i, j] + deviation_enec_fjk['wj_fjk_fij_5'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_fik_6'][i, j] + deviation_enec_fjk['wk_fjk_fik_6'][i, j]) == 0):
            ENEC_fjk[i, j] = 0
        else: ENEC_fjk[i, j] = 1
print(ENEC_fjk)

situation_initiale = 'fjk'
graph_matrice_e(structure,n_b_ci, n_b_cj, titre_b, situation_initiale,deviation_enec_fjk)
graph_matrice(structure,n_b_ci, n_b_cj, ENEC_fjk, "ENEC_$f_{jk}$")








##################################### Déviation de coalition - multilatérale #########################################
structure = "multilatérale"
# Déviation coalition - F
vect_dev_m_F = [w_quo, w_quo, w_quo, w_mij, w_mik, w_mjk]
symbole_dev = ['quo', 'quo', 'quo', 'mij', 'mik', 'mjk']
vect_pays = [['i','j'], ['i','k'], ['j','k'],['i', 'j'],['i','k'], ['j','k']]
vect_cre_m_F = [[[deviation_mik['wi_mik_quo']],[deviation_mjk['wj_mjk_quo']]],
              [[deviation_mij['wi_mij_quo']],[deviation_mjk['wk_mjk_quo']]],
              [[deviation_mij['wj_mij_quo']],[deviation_mik['wk_mik_quo']]],
              [[dev_inverse_mij['wi_quo_mij']], [dev_inverse_mij['wj_quo_mij']]],
              [[dev_inverse_mik['wi_quo_mik']], [dev_inverse_mik['wk_quo_mik']]],
              [[dev_inverse_mjk['wj_quo_mjk']], [dev_inverse_mjk['wk_quo_mjk']]],
              ]

ENEC_m_F = np.zeros((n_m_ci, n_m_cj))
deviation_enec_m_F = deviation_coal_matrice("mF",symbole_dev, w_F, vect_dev_m_F, vect_cre_m_F, vect_pays,
                                'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
for i in range(n_m_ci):
    for j in range(n_m_cj):
        if EN_m_F[i, j] == 0:
            ENEC_m_F[i, j] = 0
        elif ((deviation_enec_m_F['wi_mF_quo_1'][i, j] + deviation_enec_m_F['wj_mF_quo_1'][i, j]) == 0
              or (deviation_enec_m_F['wi_mF_quo_2'][i, j] + deviation_enec_m_F['wk_mF_quo_2'][i, j]) == 0
              or (deviation_enec_m_F['wj_mF_quo_3'][i, j] + deviation_enec_m_F['wk_mF_quo_3'][i, j]) == 0
              or (deviation_enec_m_F['wi_mF_mij_4'][i, j] + deviation_enec_m_F['wj_mF_mij_4'][i, j]) == 0
              or (deviation_enec_m_F['wi_mF_mik_5'][i, j] + deviation_enec_m_F['wk_mF_mik_5'][i, j]) == 0
              or (deviation_enec_m_F['wj_mF_mjk_6'][i, j] + deviation_enec_m_F['wk_mF_mjk_6'][i, j]) == 0):
            ENEC_m_F[i, j] = 0
        else: ENEC_m_F[i, j] = 1
print(ENEC_m_F)
inverse_ENEC_m_F = 1 - ENEC_m_F

situation_initiale = 'F'
graph_matrice_e(structure,n_m_ci, n_m_cj, titre_m, situation_initiale,deviation_enec_m_F)
graph_matrice(structure,n_m_ci, n_m_cj, ENEC_m_F, "ENEC_F")


# Déviation coalition - statu quo
vect_dev_m_quo = [w_mij, w_mik, w_mjk, w_F]
symbole_dev = ['mij', 'mik', 'mjk', 'mF']
vect_pays = [['i','j'], ['i','k'], ['j','k'], ['i', 'j', 'k']]
vect_cre_m_quo = [[[dev_inverse_mij['wi_quo_mij']],[dev_inverse_mij['wj_quo_mij']]],
                [[dev_inverse_mik['wi_quo_mik']],[dev_inverse_mik['wk_quo_mik']]],
                [[dev_inverse_mjk['wj_quo_mjk']],[dev_inverse_mjk['wk_quo_mjk']]],
                [[inverse_ENEC_m_F], [inverse_ENEC_m_F],[inverse_ENEC_m_F]]
                ]

ENEC_m_quo = np.zeros((n_m_ci, n_m_cj))
deviation_enec_m_quo = deviation_coal_matrice('quo', symbole_dev, w_quo, vect_dev_m_quo, vect_cre_m_quo,
                            vect_pays, 'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
for i in range(n_m_ci):
    for j in range(n_m_cj):
        if ((deviation_enec_m_quo['wi_quo_mij_1'][i, j] + deviation_enec_m_quo['wj_quo_mij_1'][i, j]) == 0
              or (deviation_enec_m_quo['wi_quo_mik_2'][i, j] + deviation_enec_m_quo['wk_quo_mik_2'][i, j]) == 0
              or (deviation_enec_m_quo['wj_quo_mjk_3'][i, j] + deviation_enec_m_quo['wk_quo_mjk_3'][i, j]) == 0
              or (deviation_enec_m_quo['wj_quo_mF_4'][i, j] + deviation_enec_m_quo['wk_quo_mF_4'][i, j] +
                  deviation_enec_m_quo['wj_quo_mF_4'][i, j]) == 0):
            ENEC_m_quo[i, j] = 0
        else: ENEC_m_quo[i, j] = 1
print(ENEC_m_quo)

situation_initiale = 'statu quo'
graph_matrice_e(structure,n_m_ci, n_m_cj, titre_m, situation_initiale,deviation_enec_m_quo)
graph_matrice(structure,n_m_ci, n_m_cj, ENEC_m_quo, "ENEC_$\phi$")



# Déviation coalition - mij
dev_manquant_mij, dev_manquant_inverse_mij = deviation_uni_matrice("mij",['mik', 'mjk'],
                                                    w_mij,[w_mik, w_mjk],  ['i','j','k'],
                                            'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
dev_manquant_mik, dev_manquant_inverse_mik = deviation_uni_matrice("mik",['mij', 'mjk'],
                                                w_mik,[w_mij, w_mjk],  ['i','j','k'],
                                            'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)

vect_dev_mij = [w_F, w_mik, w_mjk]
symbole_dev = ['mF', 'mik', 'mjk']
vect_pays = [['i', 'j', 'k'], ['i','k'], ['j','k']]
vect_cre_mij = [[[inverse_ENEC_m_F],[inverse_ENEC_m_F], [inverse_ENEC_m_F]],
                [[dev_manquant_mij['wi_mij_mik'], dev_inverse_mik['wi_quo_mik']], [dev_inverse_mik['wk_quo_mik']]],
                [[dev_manquant_mij['wj_mij_mjk'], dev_inverse_mjk['wj_quo_mjk']], [dev_inverse_mjk['wk_quo_mjk']]],
                ]

ENEC_mij = np.zeros((n_m_ci, n_m_cj))
deviation_enec_mij = deviation_coal_matrice('mij', symbole_dev, w_mij, vect_dev_mij, vect_cre_mij,
                                vect_pays, 'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
for i in range(n_m_ci):
    for j in range(n_m_cj):
        if EN_mij[i, j] == 0:
            ENEC_mij[i, j] = 0
        elif ((deviation_enec_mij['wi_mij_mF_1'][i, j] + deviation_enec_mij['wj_mij_mF_1'][i, j] +
               deviation_enec_mij['wk_mij_mF_1'][i, j]) == 0
              or (deviation_enec_mij['wi_mij_mik_2'][i, j] + deviation_enec_mij['wk_mij_mik_2'][i, j]) == 0
              or (deviation_enec_mij['wj_mij_mjk_3'][i, j] + deviation_enec_mij['wk_mij_mjk_3'][i, j]) == 0):
            ENEC_mij[i, j] = 0
        else: ENEC_mij[i, j] = 1
print(ENEC_mij)

situation_initiale = 'mij'
graph_matrice_e(structure,n_m_ci, n_m_cj, titre_m, situation_initiale,deviation_enec_mij)
graph_matrice(structure,n_m_ci, n_m_cj, ENEC_mij, "ENEC_$m_{ij}$")


# Déviation coalition - mik
vect_dev_mik = [w_F, w_mij, w_mjk]
symbole_dev = ['mF', 'mij', 'mjk']
vect_pays = [['i', 'j', 'k'], ['i','j'], ['j','k']]
vect_cre_mik = [[[inverse_ENEC_m_F],[inverse_ENEC_m_F], [inverse_ENEC_m_F]],
                [[dev_manquant_mik['wi_mik_mij'], dev_inverse_mij['wi_quo_mij']], [dev_inverse_mij['wj_quo_mij']]],
                [[dev_inverse_mjk['wj_quo_mjk']], [dev_manquant_mik['wk_mik_mjk'], dev_inverse_mjk['wk_quo_mjk']]],
                ]

ENEC_mik = np.zeros((n_m_ci, n_m_cj))
deviation_enec_mik = deviation_coal_matrice('mik', symbole_dev, w_mik, vect_dev_mik, vect_cre_mik,
                            vect_pays, 'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
for i in range(n_m_ci):
    for j in range(n_m_cj):
        if EN_mik[i, j] == 0:
            ENEC_mik[i, j] = 0
        elif ((deviation_enec_mik['wi_mik_mF_1'][i, j] + deviation_enec_mik['wj_mik_mF_1'][i, j] +
               deviation_enec_mik['wk_mik_mF_1'][i, j]) == 0
              or (deviation_enec_mik['wi_mik_mij_2'][i, j] + deviation_enec_mik['wj_mik_mij_2'][i, j]) == 0
              or (deviation_enec_mik['wj_mik_mjk_3'][i, j] + deviation_enec_mik['wk_mik_mjk_3'][i, j]) == 0):
            ENEC_mik[i, j] = 0
        else: ENEC_mik[i, j] = 1
print(ENEC_mik)

situation_initiale = 'mik'
graph_matrice_e(structure,n_m_ci, n_m_cj, titre_m, situation_initiale,deviation_enec_mik)
graph_matrice(structure,n_m_ci, n_m_cj, ENEC_mik, "ENEC_$m_{ik}$")



# Déviation coalition - mjk
vect_dev_mjk = [w_F, w_mij, w_mik]
symbole_dev = ['mF', 'mij', 'mik']
vect_pays = [['i', 'j', 'k'], ['i','j'], ['i','k']]
vect_cre_mjk = [[[inverse_ENEC_m_F],[inverse_ENEC_m_F], [inverse_ENEC_m_F]],
                [[dev_inverse_mij['wi_quo_mij']], [dev_manquant_inverse_mij['wj_mjk_mij'],
                                                   dev_inverse_mij['wj_quo_mij']]],
                [[dev_inverse_mik['wi_quo_mik']], [dev_manquant_inverse_mik['wk_mjk_mik'],
                                                   dev_inverse_mik['wk_quo_mik']]],
                ]

ENEC_mjk = np.zeros((n_m_ci, n_m_cj))
deviation_enec_mjk = deviation_coal_matrice('mjk', symbole_dev, w_mjk, vect_dev_mjk, vect_cre_mjk,
                        vect_pays, 'multilatérale', n_m_ci, n_m_cj, c_values_m_ci, c_values_m_cj)
for i in range(n_m_ci):
    for j in range(n_m_cj):
        if EN_mjk[i, j] == 0:
            ENEC_mjk[i, j] = 0
        elif ((deviation_enec_mjk['wi_mjk_mF_1'][i, j] + deviation_enec_mjk['wj_mjk_mF_1'][i, j] +
               deviation_enec_mjk['wk_mjk_mF_1'][i, j]) == 0
              or (deviation_enec_mjk['wi_mjk_mij_2'][i, j] + deviation_enec_mjk['wj_mjk_mij_2'][i, j]) == 0
              or (deviation_enec_mjk['wi_mjk_mik_3'][i, j] + deviation_enec_mjk['wk_mjk_mik_3'][i, j]) == 0):
            ENEC_mjk[i, j] = 0
        else: ENEC_mjk[i, j] = 1
print(ENEC_mjk)

situation_initiale = 'mjk'
graph_matrice_e(structure,n_m_ci, n_m_cj, titre_m, situation_initiale,deviation_enec_mjk)
graph_matrice(structure,n_m_ci, n_m_cj, ENEC_mjk, "ENEC_$m_{jk}$")


matrice_presentation = ENEC_F_final + ENEC_m_F
print(matrice_presentation)
graph_matrice_presentation(n_m_ci, n_m_cj, matrice_presentation, "")