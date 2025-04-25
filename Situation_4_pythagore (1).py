import sys

from sympy.conftest import sp

sys.path.append("/fonctions_projet.py")
sys.path.append("/fonctions_mj.py")

from fonctions_projet import *

############################# situation 4 pythagore #################################
x = 0.1

p_values = np.arange(0, 1+x, x)
n_p = len(p_values)

# Tarifs optimaux
s_t_o_b_4_vect, s_t_o_b_4_dict    = tarifs_optimaux_bilateraux(1, 1, 1,c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))

s_t_o_m_4_vect, s_t_o_m_4_dict    = tarifs_optimaux_multilateraux(1, 1, 1,c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))

# Prix d'equilibre
p_equ_vect, p_equ_dict = prix_equilibre(a, 1, 1, 1, s_t_o_b_4_dict["ti_quo"], s_t_o_b_4_dict["tj_quo"],
                                        s_t_o_b_4_dict["tk_quo"], s_t_o_b_4_dict["ti_fij"], s_t_o_b_4_dict["tj_fij"]
                                        , s_t_o_b_4_dict["ti_fik"], s_t_o_b_4_dict["tk_fik"], s_t_o_b_4_dict["tj_fjk"],
                                        s_t_o_b_4_dict["tk_fjk"], s_t_o_m_4_dict["ti_mij"]
                                        , s_t_o_m_4_dict["tj_mij"], s_t_o_m_4_dict["ti_mik"], s_t_o_m_4_dict["tk_mik"],
                                        s_t_o_m_4_dict["tj_mjk"], s_t_o_m_4_dict["tk_mjk"]
                                        , s_t_o_b_4_dict["tj_ih"], s_t_o_b_4_dict["tk_ih"], s_t_o_b_4_dict["ti_jh"],
                                        s_t_o_b_4_dict["tk_jh"], s_t_o_b_4_dict["ti_kh"], s_t_o_b_4_dict["tj_kh"]
                                        ,c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))




# Exportations
x_quo = export_i(a, 1, 1, 1, p_equ_dict["pi_i_quo"], p_equ_dict["pj_j_quo"], p_equ_dict["pk_k_quo"], s_t_o_b_4_dict["ti_quo"], s_t_o_b_4_dict["ti_quo"], s_t_o_b_4_dict["tj_quo"], s_t_o_b_4_dict["tj_quo"],s_t_o_b_4_dict["tk_quo"], s_t_o_b_4_dict["tk_quo"],c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
x_fij = export_i(a, 1, 1, 1, p_equ_dict["pi_i_fij"], p_equ_dict["pj_j_fij"], p_equ_dict["pk_k_fij"], 0,s_t_o_b_4_dict["ti_fij"], 0, s_t_o_b_4_dict["tj_fij"], s_t_o_b_4_dict["tk_quo"],s_t_o_b_4_dict["tk_quo"],c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
x_fik = export_i(a, 1, 1, 1, p_equ_dict["pi_i_fik"], p_equ_dict["pj_j_fik"], p_equ_dict["pk_k_fik"], s_t_o_b_4_dict["ti_fik"], 0, s_t_o_b_4_dict["tj_quo"], s_t_o_b_4_dict["tj_quo"], 0,s_t_o_b_4_dict["tk_fik"],c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
x_fjk = export_i(a, 1, 1, 1, p_equ_dict["pi_i_fjk"], p_equ_dict["pj_j_fjk"], p_equ_dict["pk_k_fjk"], s_t_o_b_4_dict["ti_quo"], s_t_o_b_4_dict["ti_quo"], s_t_o_b_4_dict["tj_fjk"], 0, s_t_o_b_4_dict["tk_fjk"], 0,c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))


print("")
x_mij = export_i(a, 1, 1, 1, p_equ_dict["pi_i_mij"], p_equ_dict["pj_j_mij"], p_equ_dict["pk_k_mij"], s_t_o_m_4_dict["ti_mij"], s_t_o_m_4_dict["ti_mij"], s_t_o_m_4_dict["tj_mij"], s_t_o_m_4_dict["tj_mij"], s_t_o_m_4_dict["tk_quo"], s_t_o_m_4_dict["tk_quo"],c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
x_mik = export_i(a, 1, 1, 1, p_equ_dict["pi_i_mik"], p_equ_dict["pj_j_mik"], p_equ_dict["pk_k_mik"], s_t_o_m_4_dict["ti_mik"], s_t_o_m_4_dict["ti_mik"], s_t_o_m_4_dict["tj_quo"], s_t_o_m_4_dict["tj_quo"], s_t_o_m_4_dict["tk_mik"], s_t_o_m_4_dict["tk_mik"],c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
x_mjk = export_i(a, 1, 1, 1, p_equ_dict["pi_i_mjk"], p_equ_dict["pj_j_mjk"], p_equ_dict["pk_k_mjk"], s_t_o_m_4_dict["ti_quo"], s_t_o_m_4_dict["ti_quo"], s_t_o_m_4_dict["tj_mjk"], s_t_o_m_4_dict["tj_mjk"],
                 s_t_o_m_4_dict["tk_mjk"], s_t_o_m_4_dict["tk_mjk"],c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
print("")
x_ih = export_i(a, 1, 1, 1, p_equ_dict["pi_i_ih"], p_equ_dict["pj_j_ih"], p_equ_dict["pk_k_ih"], 0, 0, 0, s_t_o_b_4_dict["tj_ih"], 0, s_t_o_b_4_dict["tk_ih"],c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
x_jh = export_i(a, 1, 1, 1, p_equ_dict["pi_i_jh"], p_equ_dict["pj_j_jh"], p_equ_dict["pk_k_jh"], 0, s_t_o_b_4_dict["ti_jh"], 0, 0, s_t_o_b_4_dict["tk_jh"], 0,c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
x_kh = export_i(a, 1, 1, 1, p_equ_dict["pi_i_kh"], p_equ_dict["pj_j_kh"], p_equ_dict["pk_k_kh"], s_t_o_b_4_dict["ti_kh"], 0, s_t_o_b_4_dict["tj_kh"], 0, 0, 0,c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
x_F = export_i(a, 1, 1, 1, p_equ_dict["pi_i_F"], p_equ_dict["pj_j_F"], p_equ_dict["pk_k_F"], 0, 0, 0, 0, 0, 0,c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c))
X_b, X_m = export_final(x_quo, x_fij, x_fik, x_fjk, x_mij, x_mik, x_mjk, x_ih, x_jh, x_kh, x_F)

print("*Mettre les équations en fonction de c seulement")

s_t_o_b_4_all,X_b_all,tarif_b_all,s_t_o_m_4_all,X_m_all,tarif_m_all = mettre_en_fonction_de_c(s_t_o_b_4_dict, s_t_o_m_4_dict, s_t_o_b_4_vect, s_t_o_m_4_vect, X_b, X_m, p_values)


# Définir la variable symbolique
print("*Trier des équations et les solutionner")
variable = sp.symbols('c')


liste_solutions_bilaterales_c, liste_problematique_b,liste_eq_sans_solution_b,liste_solutions_multilaterales_c, liste_problematique_m,liste_eq_sans_solution_m = trier_les_equations_et_solver(s_t_o_b_4_all,X_b_all,tarif_b_all,s_t_o_m_4_all, X_m_all,tarif_m_all)





############################ bornes bilatérales ####################################
print("Trouver les bornes")
# Affichage du type et des éléments dans liste_solutions_bilaterales_c pour débogage
print("Vérification des inégalités dans liste_solutions_bilaterales_c :")
for sublist in liste_solutions_bilaterales_c:
    print(f"Liste d'inégalités : {sublist}, Type : {type(sublist)}")
    if isinstance(sublist, list):
        for ineq in sublist:
            print(f"Inégalité : {ineq}, Type : {type(ineq)}")
    elif isinstance(sublist, sp.And):
        print(f"Conjonction d'inégalités (And) : {sublist}, Type : {type(sublist)}")

# Normaliser les inégalités
inegalites_normalisees_bilaterales= []
for sublist in liste_solutions_bilaterales_c:
    if isinstance(sublist, list):
        inegalites_normalisees_bilaterales.extend(normaliser_inegalites(sublist))
    elif isinstance(sublist, sp.And):
        inegalites_normalisees_bilaterales.extend(normaliser_inegalites([sublist]))


bornes_b_min, bornes_b_max = extraire_bornes_bilaterales(inegalites_normalisees_bilaterales)



############################ bornes multilatérales ####################################


    
# Affichage du type et des éléments dans liste_solutions_multilaterales_c pour débogage
print("Vérification des inégalités dans liste_solutions_multilaterales_c :")
for sublist in liste_solutions_multilaterales_c:
    print(f"Liste d'inégalités : {sublist}, Type : {type(sublist)}")
    if isinstance(sublist, list):
        for ineq in sublist:
            print(f"Inégalité : {ineq}, Type : {type(ineq)}")
    elif isinstance(sublist, sp.And):
        print(f"Conjonction d'inégalités (And) : {sublist}, Type : {type(sublist)}")

# Normaliser les inégalités
inegalites_normalisees_multilaterales= []
for sublist in liste_solutions_multilaterales_c:
    if isinstance(sublist, list):
        inegalites_normalisees_multilaterales.extend(normaliser_inegalites(sublist))
    elif isinstance(sublist, sp.And):
        inegalites_normalisees_multilaterales.extend(normaliser_inegalites([sublist]))



bornes_m_min, bornes_m_max = extraire_bornes_multilaterales(inegalites_normalisees_multilaterales)


print("Les biens-êtres optimaux")
########## LES BIENS-ÊTRES OPTIMAUX ##########
# Bien-être optimaux - Structure bilaterale
print("##### Statu quo ####")
w_quo = w_optimaux(1, 1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), s_t_o_b_4_dict["ti_quo"], s_t_o_b_4_dict["tj_quo"], s_t_o_b_4_dict["ti_quo"], s_t_o_b_4_dict["tk_quo"], s_t_o_b_4_dict["tj_quo"], s_t_o_b_4_dict["tk_quo"],"quo")
print("")

print("#### fij ####")
w_fij = w_optimaux(1, 1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), 0, 0, s_t_o_b_4_dict["ti_fij"], s_t_o_b_4_dict["tk_quo"], s_t_o_b_4_dict["tj_fij"], s_t_o_b_4_dict["tk_quo"], "fij")
print("")

print("##### fik ####")
w_fik = w_optimaux(1,1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), s_t_o_b_4_dict["ti_fik"], s_t_o_b_4_dict["tj_quo"], 0, 0, s_t_o_b_4_dict["tj_quo"], s_t_o_b_4_dict["tk_fik"], "fik")
print("")

print("##### fjk ####")
w_fjk = w_optimaux(1,1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), s_t_o_b_4_dict["ti_quo"], s_t_o_b_4_dict["tj_fjk"], s_t_o_b_4_dict["ti_quo"], s_t_o_b_4_dict["tk_fjk"], 0, 0, "fjk")
print("")

print("##### ih ####")
w_ih = w_optimaux(1,1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), 0, 0, 0, 0, s_t_o_b_4_dict["tj_ih"], s_t_o_b_4_dict["tk_ih"], "ih")
print("")

print("##### jh ####")
w_jh = w_optimaux(1,1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), 0, 0, s_t_o_b_4_dict["ti_jh"], s_t_o_b_4_dict["tk_jh"], 0, 0, "jh")
print("")

print("##### kh ####")
w_kh = w_optimaux(1,1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), s_t_o_b_4_dict["ti_kh"], s_t_o_b_4_dict["tj_kh"], 0, 0, 0, 0, "kh")
print("")

print("##### F ####")
w_F = w_optimaux(1,1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), 0, 0, 0, 0, 0, 0, "F")
print("")

# Bien-être optimaux - Structure multilaterale
print("##### mij ####")
w_mij = w_optimaux(1,1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), s_t_o_m_4_dict["ti_mij"], s_t_o_m_4_dict["tj_mij"], s_t_o_m_4_dict["ti_mij"], s_t_o_m_4_dict["tk_quo"], s_t_o_m_4_dict["tj_mij"], s_t_o_m_4_dict["tk_quo"], "mij")
print("")

print("##### mik ####")
w_mik = w_optimaux(1,1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), s_t_o_m_4_dict["ti_mik"], s_t_o_m_4_dict["tj_quo"], s_t_o_m_4_dict["ti_mik"], s_t_o_m_4_dict["tk_mik"], s_t_o_m_4_dict["tj_quo"], s_t_o_m_4_dict["tk_mik"], "mik")
print("")

print("##### mjk ####")
w_mjk = w_optimaux(1,1, 1, 1, c, c, sqrt((p*c)**2+c**2), sqrt((p*c)**2+c**2), (p*c), (p*c), s_t_o_m_4_dict["ti_quo"], s_t_o_m_4_dict["tj_mjk"], s_t_o_m_4_dict["ti_quo"], s_t_o_m_4_dict["tk_mjk"], s_t_o_m_4_dict["tj_mjk"], s_t_o_m_4_dict["tk_mjk"], "mjk")
print("")



######################################### Déviation unilatérale - Bilatérale #######################################################
print("les deviations unilaterales")
print(bornes_b_min)
print(bornes_b_max)
print(p_values)
print(n_p)


bij = float(bornes_b_max.evalf())
bjk = float(max(p_values) * bij)

print(bornes_m_min)
print(bornes_m_max)
print(p_values)
print(n_p)


mij = float(bornes_m_max.evalf())
mjk = float(max(p_values) * mij)


c_values_b_cij = np.arange(0, round(mij+x,2), x)
c_values_b_cjk = np.arange(0, round(mjk+x,2), x)
n_b_cij = len(c_values_b_cij)
n_b_cjk = len(c_values_b_cjk)



c_values_m_cij = np.arange(0, round(mij+x,2), x)
c_values_m_cjk = np.arange(0, round(mjk+x,2), x)
n_m_cij = len(c_values_m_cij)
n_m_cjk = len(c_values_m_cjk)

vect_pays = ['i','j','k']
titre_b = 'unilatérale bilatérale'
titre_m = 'unilatérale multilatérale'


# Déviation unilatérale - fij
vect_dev_fij = [w_quo]
symbole_dev = ['quo']
EN_fij = np.zeros((n_b_cij, n_p))
deviation_fij, dev_inverse_fij = deviation_uni_matrice_pythagore("fij",symbole_dev, w_fij, vect_dev_fij,vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if deviation_fij['wi_fij_quo'][i,j] == 2 and deviation_fij['wj_fij_quo'][i,j] == 2:
            EN_fij[i,j] = 2
        elif deviation_fij['wi_fij_quo'][i,j] == 1 and deviation_fij['wj_fij_quo'][i,j] == 1:
            EN_fij[i,j] = 1
        else: EN_fij[i,j] = 0
print(EN_fij)
inverse_fij = 1- EN_fij

situation_initiale = 'fij'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale, deviation_fij)
graph_matrice_pythagore(n_b_cij, n_p, EN_fij, "EN_fij")





# Déviation unilatérale - fik
vect_dev_fik = [w_quo]
symbole_dev = ['quo']
EN_fik = np.zeros((n_b_cij, n_p))
deviation_fik, dev_inverse_fik = deviation_uni_matrice_pythagore("fik",symbole_dev, w_fik, vect_dev_fik,vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if deviation_fik['wi_fik_quo'][i,j] == 2 and deviation_fik['wk_fik_quo'][i,j] == 2:
            EN_fik[i,j] = 2
        elif deviation_fik['wi_fik_quo'][i,j] == 1 and deviation_fik['wk_fik_quo'][i,j] == 1:
            EN_fik[i,j] = 1
        else: EN_fik[i,j] = 0
print(EN_fik)
inverse_fik = 1 - EN_fik

situation_initiale = 'fik'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_fik)
graph_matrice_pythagore(n_b_cij, n_p, EN_fik, "EN_fik")



# Déviation unilatérale - fjk
vect_dev_fjk = [w_quo]
symbole_dev = ['quo']
EN_fjk = np.zeros((n_b_cij, n_p))
deviation_fjk, dev_inverse_fjk = deviation_uni_matrice_pythagore("fjk",symbole_dev, w_fjk, vect_dev_fjk, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if deviation_fjk['wj_fjk_quo'][i,j] == 2 and deviation_fjk['wk_fjk_quo'][i,j] == 2:
            EN_fjk[i,j] = 2
        elif deviation_fjk['wj_fjk_quo'][i,j] == 1 and deviation_fjk['wk_fjk_quo'][i,j] == 1:
            EN_fjk[i,j] = 1
        else: EN_fjk[i,j] = 0
print(EN_fjk)
inverse_fjk = 1 - EN_fjk

situation_initiale = 'fjk'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_fjk)
graph_matrice_pythagore(n_b_cij, n_p, EN_fjk, "EN_fjk")



# Déviation unilatérale - ih
vect_dev_ih = [w_fij, w_fik, w_quo]
symbole_dev = ['fij', 'fik', 'quo']
EN_ih = np.zeros((n_b_cij, n_p))
deviation_ih, dev_inverse_ih = deviation_uni_matrice_pythagore("ih",symbole_dev, w_ih, vect_dev_ih, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if (deviation_ih['wi_ih_fij'][i,j] == 2 and deviation_ih['wi_ih_fik'][i,j] == 2 and deviation_ih['wi_ih_quo'][i,j] == 2
                and deviation_ih['wj_ih_fik'][i,j] == 2 and deviation_ih['wk_ih_fij'][i,j] == 2):
            EN_ih[i,j] = 2
        elif (deviation_ih['wi_ih_fij'][i,j] == 1 and deviation_ih['wi_ih_fik'][i,j] == 1 and deviation_ih['wi_ih_quo'][i,j] == 1
                and deviation_ih['wj_ih_fik'][i,j] == 1 and deviation_ih['wk_ih_fij'][i,j] == 1):
            EN_ih[i,j] = 1
        else: EN_ih[i,j] = 0
print(EN_ih)
inverse_ih = 1 - EN_ih

situation_initiale = 'ih'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_ih)
graph_matrice_pythagore(n_b_cij, n_p, EN_ih, "EN_ih")


# Déviation unilatérale - jh
vect_dev_jh = [w_fjk, w_fij, w_quo]
symbole_dev = ['fjk', 'fij', 'quo']
EN_jh = np.zeros((n_b_cij, n_p))
deviation_jh, dev_inverse_jh= deviation_uni_matrice_pythagore("jh",symbole_dev, w_jh, vect_dev_jh, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if (deviation_jh['wj_jh_fjk'][i,j] == 2 and deviation_jh['wj_jh_fij'][i,j] == 2 and deviation_jh['wj_jh_quo'][i,j] == 2
                and deviation_jh['wk_jh_fij'][i,j] == 2 and deviation_jh['wi_jh_fjk'][i,j] == 2):
            EN_jh[i,j] = 2
        elif (deviation_jh['wj_jh_fjk'][i,j] == 1 and deviation_jh['wj_jh_fij'][i,j] == 1 and deviation_jh['wj_jh_quo'][i,j] == 1
                and deviation_jh['wk_jh_fij'][i,j] == 1 and deviation_jh['wi_jh_fjk'][i,j] == 1):
            EN_jh[i,j] = 1
        else: EN_jh[i,j] = 0
print(EN_jh)
inverse_jh = 1 - EN_jh
situation_initiale = 'jh'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_jh)
graph_matrice_pythagore(n_b_cij, n_p, EN_jh, "EN_jh")



# Déviation unilatérale - kh
vect_dev_kh = [w_fjk, w_fik, w_quo]
symbole_dev = ['fjk', 'fik', 'quo']
EN_kh = np.zeros((n_b_cij, n_p))
deviation_kh, dev_inverse_kh = deviation_uni_matrice_pythagore("kh",symbole_dev, w_kh, vect_dev_kh, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if (deviation_kh['wk_kh_fjk'][i,j] == 2 and deviation_kh['wk_kh_fik'][i,j] == 2 and deviation_kh['wk_kh_quo'][i,j] == 2 
                and deviation_kh['wj_kh_fik'][i,j] == 2 and deviation_kh['wi_kh_fjk'][i,j] == 2):
            EN_kh[i,j] = 2            
        elif (deviation_kh['wk_kh_fjk'][i,j] == 1 and deviation_kh['wk_kh_fik'][i,j] == 1 and deviation_kh['wk_kh_quo'][i,j] == 1 
                and deviation_kh['wj_kh_fik'][i,j] == 1 and deviation_kh['wi_kh_fjk'][i,j] == 1):
            EN_kh[i,j] = 1
        else: EN_kh[i,j] = 0
print(EN_kh)
inverse_kh = 1 - EN_kh

situation_initiale = 'kh'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_kh)
graph_matrice_pythagore(n_b_cij, n_p, EN_kh, "EN_kh")


# Déviation unilatérale - F
vect_dev_F = [w_jh, w_kh, w_fjk, w_ih, w_fik, w_fij]
symbole_dev = ['jh', 'kh', 'fjk', 'ih', 'fik', 'fij']
EN_F = np.zeros((n_b_cij, n_p))
deviation_F, dev_inverse_F = deviation_uni_matrice_pythagore("F",symbole_dev, w_F, vect_dev_F, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if (deviation_F['wi_F_jh'][i, j] == 2 and deviation_F['wi_F_kh'][i, j] == 2 and deviation_F['wi_F_fjk'][i, j] == 2 and
                deviation_F['wj_F_ih'][i, j] == 2 and deviation_F['wj_F_kh'][i, j] == 2 and deviation_F['wj_F_fik'][i, j] == 2 and
                deviation_F['wk_F_ih'][i, j] == 2 and deviation_F['wk_F_jh'][i, j] == 2 and deviation_F['wk_F_fij'][i, j] == 2):
            EN_F[i, j] = 2
        elif (deviation_F['wi_F_jh'][i, j] == 1 and deviation_F['wi_F_kh'][i, j] == 1 and deviation_F['wi_F_fjk'][i, j] == 1 and
                deviation_F['wj_F_ih'][i, j] == 1 and deviation_F['wj_F_kh'][i, j] == 1 and deviation_F['wj_F_fik'][i, j] == 1 and
                deviation_F['wk_F_ih'][i, j] == 1 and deviation_F['wk_F_jh'][i, j] == 1 and deviation_F['wk_F_fij'][i, j] == 1):
            EN_F[i, j] = 1
        else: EN_F[i, j] = 0
print(EN_F)
inverse_F = 1 - EN_F

situation_initiale = 'F'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_F)
graph_matrice_pythagore(n_b_cij, n_p, EN_F, "EN_F")


######################################### Déviation unilatérale - multilatérale #######################################################



# Déviation unilatérale - mij
vect_dev_mij = [w_quo]
symbole_dev = ['quo']
EN_mij = np.zeros((n_m_cij, n_p))
deviation_mij, dev_inverse_mij = deviation_uni_matrice_pythagore("mij",symbole_dev, w_mij, vect_dev_mij,vect_pays, 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
for i in range(n_m_cij):
    for j in range(n_p):
        if deviation_mij['wi_mij_quo'][i, j] == 2 and deviation_mij['wj_mij_quo'][i, j] == 2:
            EN_mij[i, j] = 2
        elif deviation_mij['wi_mij_quo'][i, j] == 1 and deviation_mij['wj_mij_quo'][i, j] == 1:
            EN_mij[i, j] = 1
        else: EN_mij[i, j] = 0
print(EN_mij)
inverse_mij = 1 - EN_mij

situation_initiale = 'mij'
graph_matrice_e_pythagore(n_m_cij, n_p, titre_m, situation_initiale,deviation_mij)
graph_matrice_pythagore(n_m_cij, n_p, EN_mij, "EN_mij")

# Déviation unilatérale - mik
vect_dev_mik = [w_quo]
symbole_dev = ['quo']
EN_mik = np.zeros((n_m_cij, n_p))
deviation_mik, dev_inverse_mik = deviation_uni_matrice_pythagore("mik",symbole_dev, w_mik, vect_dev_mik,vect_pays, 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
for i in range(n_m_cij):
    for j in range(n_p):
        if deviation_mik['wi_mik_quo'][i, j] == 2 and deviation_mik['wk_mik_quo'][i, j] == 2:
            EN_mik[i, j] = 2
        if deviation_mik['wi_mik_quo'][i, j] == 1 and deviation_mik['wk_mik_quo'][i, j] == 1:
            EN_mik[i, j] = 1
        else: EN_mik[i, j] = 0
print(EN_mik)
inverse_mik = 1 - EN_mik

situation_initiale = 'mik'
graph_matrice_e_pythagore(n_m_cij, n_p, titre_m, situation_initiale,deviation_mik)
graph_matrice_pythagore(n_m_cij, n_p, EN_mik, "EN_mik")


# Déviation unilatérale - mjk
vect_dev_mjk = [w_quo]
symbole_dev = ['quo']
EN_mjk = np.zeros((n_m_cij, n_p))
deviation_mjk, dev_inverse_mjk = deviation_uni_matrice_pythagore("mjk",symbole_dev, w_mjk, vect_dev_mjk, vect_pays, 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
for i in range(n_m_cij):
    for j in range(n_p):
        if deviation_mjk['wj_mjk_quo'][i, j] == 2 and deviation_mjk['wk_mjk_quo'][i, j] == 2:
            EN_mjk[i, j] = 2
        elif deviation_mjk['wj_mjk_quo'][i, j] == 1 and deviation_mjk['wk_mjk_quo'][i, j] == 1:
            EN_mjk[i, j] = 1
        else: EN_mjk[i, j] = 0
print(EN_mjk)
inverse_mjk = 1 - EN_mjk

situation_initiale = 'mjk'
graph_matrice_e_pythagore(n_m_cij, n_p, titre_m, situation_initiale,deviation_mjk)
graph_matrice_pythagore(n_m_cij, n_p, EN_mjk, "EN_mjk")



# Déviation unilatérale - F
vect_dev_F = [w_mjk, w_mik, w_mij]
symbole_dev = ['mjk', 'mik', 'mij']
EN_m_F = np.zeros((n_m_cij, n_p))
deviation_m_F, dev_inverse_m_F = deviation_uni_matrice_pythagore("F",symbole_dev, w_F, vect_dev_F, vect_pays, 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
for i in range(n_m_cij):
    for j in range(n_p):
        if deviation_m_F['wi_F_mjk'][i, j] == 2 and deviation_m_F['wj_F_mik'][i, j] == 2 and deviation_m_F['wk_F_mij'][i, j] == 2:
            EN_m_F[i, j] = 2
        elif deviation_m_F['wi_F_mjk'][i, j] == 1 and deviation_m_F['wj_F_mik'][i, j] == 1 and deviation_m_F['wk_F_mij'][i, j] == 1:
            EN_m_F[i, j] = 1
        else: EN_m_F[i, j] = 0
print(EN_m_F)
inverse_m_F = 1 - EN_m_F

situation_initiale = 'F'
graph_matrice_e_pythagore(n_m_cij, n_p, titre_m, situation_initiale,deviation_m_F)
graph_matrice_pythagore(n_m_cij, n_p, EN_m_F, "EN_m_F")



##################################### Déviation de coalition - Bilatérale ##############################################
print("Les deviations de coalition")
titre_b = 'de coalition bilatérale'
titre_m = 'de coalition multilatérale'

#Déviation coalition - F
vect_dev_F = [w_quo, w_quo, w_quo, w_fij, w_fik, w_fjk]
symbole_dev = ['quo', 'quo', 'quo', 'fij', 'fik', 'fjk']
vect_pays = [['i','j'], ['i','k'], ['j','k'],['i', 'j'],['i','k'], ['j','k']]
vect_cre_F = [[[deviation_fik['wi_fik_quo']],[deviation_fjk['wj_fjk_quo']]],
              [[deviation_fij['wi_fij_quo']],[deviation_fjk['wk_fjk_quo']]],
              [[deviation_fij['wj_fij_quo']],[deviation_fik['wk_fik_quo']]],
              [[deviation_ih['wi_ih_fij'],dev_inverse_fij['wi_quo_fij']], [deviation_jh['wj_jh_fij'],dev_inverse_fij['wj_quo_fij']]],
              [[deviation_ih['wi_ih_fik'],dev_inverse_fik['wi_quo_fik']], [deviation_kh['wk_kh_fik'],dev_inverse_fik['wk_quo_fik']]],
              [[deviation_jh['wj_jh_fjk'],dev_inverse_fjk['wj_quo_fjk']], [deviation_kh['wk_kh_fjk'],dev_inverse_fjk['wk_quo_fjk']]],
              ]

ENEC_F = np.zeros((n_b_cij, n_p))
deviation_enec_F = deviation_coal_matrice_pythagore("F",symbole_dev, w_F, vect_dev_F, vect_cre_F, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if EN_F[i,j] == 0:
            ENEC_F[i,j] = 0
        elif EN_F[i,j] == 2:
            ENEC_F[i,j] = 2
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
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_enec_F)
graph_matrice_pythagore(n_b_cij, n_p, ENEC_F, "ENEC_F")




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


ENEC_quo = np.zeros((n_b_cij, n_p))
deviation_enec_quo = deviation_coal_matrice_pythagore('quo', symbole_dev, w_quo, vect_dev_quo, vect_cre_quo, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if deviation_enec_quo['wi_quo_fij_1'][i,j] == 2:
            ENEC_quo[i,j] = 2
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
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_enec_quo)
graph_matrice_pythagore(n_b_cij, n_p, ENEC_quo, "ENEC_statu_quo")


# Déviation coalition - ih
# déviation manquante
dev_manquant_fij, dev_manquant_inverse_fij = deviation_uni_matrice_pythagore("fij",['fik', 'fjk'], w_fij,[w_fik, w_fjk],  ['i','j','k'], 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
dev_manquant_fik, dev_manquant_inverse_fik = deviation_uni_matrice_pythagore("fik",['fij', 'fjk'], w_fik,[w_fij, w_fjk],  ['i','j','k'],'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)

vect_dev_ih = [w_F, w_quo, w_jh, w_kh, w_fjk]
symbole_dev = ['F', 'quo', 'jh', 'kh', 'fjk']
vect_pays = [['j','k'], ['j','k'], ['j','k'],['j','k'],['j','k']]
vect_cre_ih = [[[dev_inverse_F['wj_ih_F'], dev_inverse_F['wj_kh_F'], dev_inverse_F['wj_fik_F']],[dev_inverse_F['wk_ih_F'], dev_inverse_F['wk_jh_F'], dev_inverse_F['wk_fij_F']]],
                [[deviation_fij['wj_fij_quo']],[deviation_fik['wk_fik_quo']]],
                [[dev_inverse_jh['wj_fij_jh'], dev_inverse_jh['wj_fjk_jh'], dev_inverse_jh['wj_quo_jh']],[dev_inverse_jh['wk_fij_jh'], deviation_F['wk_F_jh']]],
                [[dev_inverse_kh['wj_fik_kh'], deviation_F['wj_F_kh']], [dev_inverse_kh['wk_fik_kh'], dev_inverse_kh['wk_fjk_kh'], dev_inverse_kh['wk_quo_kh']]],
                [[deviation_jh['wj_jh_fjk'], dev_manquant_fij['wj_fij_fjk'], dev_inverse_fjk['wj_quo_fjk']],[deviation_kh['wk_kh_fjk'], dev_manquant_fik['wk_fik_fjk'], dev_inverse_fjk['wk_quo_fjk']]]
                ]

ENEC_ih = np.zeros((n_b_cij, n_p))
deviation_enec_ih = deviation_coal_matrice_pythagore('ih', symbole_dev, w_ih, vect_dev_ih, vect_cre_ih, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if EN_ih[i, j] == 2:
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
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_enec_ih)
graph_matrice_pythagore(n_b_cij, n_p, ENEC_ih, "ENEC_ih")




# Déviation coalition - jh
vect_dev_jh = [w_F, w_quo, w_ih, w_kh, w_fik]
symbole_dev = ['F', 'quo', 'ih', 'kh', 'fik']
vect_pays = [['i','k'], ['i','k'], ['i','k'],['i','k'],['i','k']]
vect_cre_jh = [[[dev_inverse_F['wi_jh_F'], dev_inverse_F['wi_kh_F'], dev_inverse_F['wi_fjk_F']],[dev_inverse_F['wk_ih_F'], dev_inverse_F['wk_jh_F'], dev_inverse_F['wk_fij_F']]],
                [[deviation_fij['wi_fij_quo']],[deviation_fjk['wk_fjk_quo']]],
                [[dev_inverse_ih['wi_fij_ih'], dev_inverse_ih['wi_fik_ih'], dev_inverse_ih['wi_quo_ih']],[dev_inverse_ih['wk_fij_ih'], deviation_F['wk_F_ih']]],
                [[dev_inverse_kh['wi_fjk_kh'], deviation_F['wi_F_kh']], [dev_inverse_kh['wk_fik_kh'], dev_inverse_kh['wk_fjk_kh'], dev_inverse_kh['wk_quo_kh']]],
                [[deviation_ih['wi_ih_fik'], dev_manquant_fij['wi_fij_fik'], dev_inverse_fik['wi_quo_fik']],[deviation_kh['wk_kh_fik'], dev_manquant_inverse_fik['wk_fjk_fik'], dev_inverse_fik['wk_quo_fik']]]
                ]

ENEC_jh = np.zeros((n_b_cij, n_p))
deviation_enec_jh = deviation_coal_matrice_pythagore('jh', symbole_dev, w_jh, vect_dev_jh, vect_cre_jh, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if EN_jh[i, j] == 2:
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
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_enec_jh)
graph_matrice_pythagore(n_b_cij, n_p, ENEC_jh, "ENEC_jh")



# Déviation coalition - kh
vect_dev_kh = [w_F, w_quo, w_ih, w_jh, w_fij]
symbole_dev = ['F', 'quo', 'ih', 'jh', 'fij']
vect_pays = [['i','j'], ['i','j'], ['i','j'],['i','j'],['i','j']]
vect_cre_kh = [[[dev_inverse_F['wi_kh_F'], dev_inverse_F['wi_jh_F'], dev_inverse_F['wi_fjk_F']],[dev_inverse_F['wj_ih_F'], dev_inverse_F['wj_kh_F'], dev_inverse_F['wj_fik_F']]],
                [[deviation_fik['wi_fik_quo']],[deviation_fjk['wj_fjk_quo']]],
                [[dev_inverse_ih['wi_fij_ih'], dev_inverse_ih['wi_fik_ih'], dev_inverse_ih['wi_quo_ih']],[dev_inverse_ih['wj_fik_ih'], deviation_F['wj_F_ih']]],
                [[dev_inverse_jh['wi_fjk_jh'], deviation_F['wi_F_jh']], [dev_inverse_jh['wj_fij_jh'], dev_inverse_jh['wj_fjk_jh'], dev_inverse_jh['wj_quo_jh']]],
                [[deviation_ih['wi_ih_fij'], dev_manquant_fik['wi_fik_fij'], dev_inverse_fij['wi_quo_fij']],[deviation_jh['wj_jh_fij'], dev_manquant_inverse_fij['wj_fjk_fij'], dev_inverse_fij['wj_quo_fij']]]
                ]

ENEC_kh = np.zeros((n_b_cij, n_p))
deviation_enec_kh = deviation_coal_matrice_pythagore('kh', symbole_dev, w_kh, vect_dev_kh, vect_cre_kh, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if EN_kh[i, j] == 2:
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
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_enec_kh)
graph_matrice_pythagore(n_b_cij, n_p, ENEC_kh, "ENEC_kh")



# Déviation coalition - fij
vect_dev_fij = [w_F, w_ih, w_jh, w_kh, w_fik, w_fjk]
symbole_dev = ['F', 'ih', 'jh', 'kh', 'fik', 'fjk']
vect_pays = [['i', 'j', 'k'], ['i','k'], ['j','k'],['i', 'j', 'k'],['i','k'], ['j','k']]
vect_cre_fij = [[[inverse_ENEC_F],[inverse_ENEC_F], [inverse_ENEC_F]],
                [[dev_inverse_ih['wi_fij_ih'], dev_inverse_ih['wi_fik_ih'], dev_inverse_ih['wi_quo_ih']],[dev_inverse_ih['wk_fij_ih']]],
                [[dev_inverse_jh['wj_fij_jh'], dev_inverse_jh['wj_fjk_jh'], dev_inverse_jh['wj_quo_jh']],[dev_inverse_jh['wk_fij_jh']]],
                [[inverse_ENEC_kh],[inverse_ENEC_kh], [inverse_ENEC_kh]],
                [[deviation_ih['wi_ih_fik'], dev_manquant_fij['wi_fij_fik'], dev_inverse_fik['wi_quo_fik']], [dev_inverse_fik['wk_quo_fik']]],
                [[deviation_jh['wj_jh_fjk'], dev_manquant_fij['wj_fij_fjk'], dev_inverse_fjk['wj_quo_fjk']], [dev_inverse_fjk['wk_quo_fjk']]],
                ]

ENEC_fij = np.zeros((n_b_cij, n_p))
deviation_enec_fij = deviation_coal_matrice_pythagore('fij', symbole_dev, w_fij, vect_dev_fij, vect_cre_fij, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if EN_fij[i, j] == 2:
            ENEC_fij[i, j] = 2
        elif EN_fij[i, j] == 0:
            ENEC_fij[i, j] = 0
        elif ((deviation_enec_fij['wi_fij_F_1'][i, j] + deviation_enec_fij['wj_fij_F_1'][i, j] + deviation_enec_fij['wk_fij_F_1'][i, j]) == 0
              or (deviation_enec_fij['wi_fij_ih_2'][i, j] + deviation_enec_fij['wk_fij_ih_2'][i, j]) == 0
              or (deviation_enec_fij['wj_fij_jh_3'][i, j] + deviation_enec_fij['wk_fij_jh_3'][i, j]) == 0
              or (deviation_enec_fij['wi_fij_kh_4'][i, j] + deviation_enec_fij['wj_fij_kh_4'][i, j] + deviation_enec_fij['wk_fij_kh_4'][i, j]) == 0
              or (deviation_enec_fij['wi_fij_fik_5'][i, j] + deviation_enec_fij['wk_fij_fik_5'][i, j]) == 0
              or (deviation_enec_fij['wj_fij_fjk_6'][i, j] + deviation_enec_fij['wk_fij_fjk_6'][i, j]) == 0):
            ENEC_fij[i, j] = 0
        else: ENEC_fij[i, j] = 1
print(ENEC_fij)

situation_initiale = 'fij'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_enec_fij)
graph_matrice_pythagore(n_b_cij, n_p, ENEC_fij, "ENEC_fij")



# Déviation coalition - fik
vect_dev_fik = [w_F, w_ih, w_jh, w_kh, w_fij, w_fjk]
symbole_dev = ['F', 'ih', 'jh', 'kh', 'fij', 'fjk']
vect_pays = [['i', 'j', 'k'], ['i','j'], ['i', 'j', 'k'], ['j','k'], ['i','j'], ['j','k']]
vect_cre_fik = [[[inverse_ENEC_F],[inverse_ENEC_F], [inverse_ENEC_F]],
                [[dev_inverse_ih['wi_fij_ih'], dev_inverse_ih['wi_fik_ih'], dev_inverse_ih['wi_quo_ih']],[dev_inverse_ih['wj_fik_ih']]],
                [[inverse_ENEC_jh],[inverse_ENEC_jh], [inverse_ENEC_jh]],
                [[dev_inverse_kh['wj_fik_kh']], [dev_inverse_kh['wk_fik_kh'], dev_inverse_kh['wk_fjk_kh'], dev_inverse_kh['wk_quo_kh']]],
                [[deviation_ih['wi_ih_fij'], dev_manquant_fik['wi_fik_fij'], dev_inverse_fij['wi_quo_fij']], [dev_inverse_fij['wj_quo_fij']]],
                [[dev_inverse_fjk['wj_quo_fjk']], [deviation_kh['wk_kh_fjk'], dev_manquant_fik['wk_fik_fjk'], dev_inverse_fjk['wk_quo_fjk']]],
                ]

ENEC_fik = np.zeros((n_b_cij, n_p))
deviation_enec_fik = deviation_coal_matrice_pythagore('fik', symbole_dev, w_fik, vect_dev_fik, vect_cre_fik, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if EN_fik[i, j] == 2:
            ENEC_fik[i, j] = 2
        elif EN_fik[i, j] == 0:
            ENEC_fik[i, j] = 0
        elif ((deviation_enec_fik['wi_fik_F_1'][i, j] + deviation_enec_fik['wj_fik_F_1'][i, j] + deviation_enec_fik['wk_fik_F_1'][i, j]) == 0
              or (deviation_enec_fik['wi_fik_ih_2'][i, j] + deviation_enec_fik['wj_fik_ih_2'][i, j]) == 0
              or (deviation_enec_fik['wi_fik_jh_3'][i, j] + deviation_enec_fik['wj_fik_jh_3'][i, j] + deviation_enec_fik['wk_fik_jh_3'][i, j]) == 0
              or (deviation_enec_fik['wj_fik_kh_4'][i, j] + deviation_enec_fik['wk_fik_kh_4'][i, j]) == 0
              or (deviation_enec_fik['wi_fik_fij_5'][i, j] + deviation_enec_fik['wj_fik_fij_5'][i, j]) == 0
              or (deviation_enec_fik['wj_fik_fjk_6'][i, j] + deviation_enec_fik['wk_fik_fjk_6'][i, j]) == 0):
            ENEC_fik[i, j] = 0
        else: ENEC_fik[i, j] = 1
print(ENEC_fik)

situation_initiale = 'fik'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_enec_fik)
graph_matrice_pythagore(n_b_cij, n_p, ENEC_fik, "ENEC_fik")


# Déviation coalition - fjk
vect_dev_fjk = [w_F, w_ih, w_jh, w_kh, w_fij, w_fik]
symbole_dev = ['F', 'ih', 'jh', 'kh', 'fij', 'fik']
vect_pays = [['i', 'j', 'k'], ['i', 'j', 'k'], ['i','j'], ['i','k'], ['i','j'], ['i','k']]
vect_cre_fjk = [[[inverse_ENEC_F],[inverse_ENEC_F], [inverse_ENEC_F]],
                [[inverse_ENEC_ih],[inverse_ENEC_ih], [inverse_ENEC_ih]],
                [[dev_inverse_jh['wi_fjk_jh']], [dev_inverse_jh['wj_fij_jh'], dev_inverse_jh['wj_fjk_jh'], dev_inverse_jh['wj_quo_jh']]],
                [[dev_inverse_kh['wi_fjk_kh']], [dev_inverse_kh['wk_fik_kh'], dev_inverse_kh['wk_fjk_kh'], dev_inverse_kh['wk_quo_kh']]],
                [[dev_inverse_fij['wi_quo_fij']], [deviation_jh['wj_jh_fij'], dev_manquant_inverse_fij['wj_fjk_fij'], dev_inverse_fij['wj_quo_fij']]],
                [[dev_inverse_fik['wi_quo_fik']], [deviation_kh['wk_kh_fik'], dev_manquant_inverse_fik['wk_fjk_fik'], dev_inverse_fik['wk_quo_fik']]],
                ]


ENEC_fjk = np.zeros((n_b_cij, n_p))
deviation_enec_fjk = deviation_coal_matrice_pythagore('fjk', symbole_dev, w_fjk, vect_dev_fjk, vect_cre_fjk, vect_pays, 'bilatérale', n_b_cij, n_p, c_values_b_cij, p_values)
for i in range(n_b_cij):
    for j in range(n_p):
        if EN_fjk[i, j] == 2:
            ENEC_fjk[i, j] = 2
        elif EN_fjk[i, j] == 0:
            ENEC_fjk[i, j] = 0
        elif ((deviation_enec_fjk['wi_fjk_F_1'][i, j] + deviation_enec_fjk['wj_fjk_F_1'][i, j] + deviation_enec_fjk['wk_fjk_F_1'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_ih_2'][i, j] + deviation_enec_fjk['wj_fjk_ih_2'][i, j] + deviation_enec_fjk['wk_fjk_ih_2'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_jh_3'][i, j] + deviation_enec_fjk['wj_fjk_jh_3'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_kh_4'][i, j] + deviation_enec_fjk['wk_fjk_kh_4'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_fij_5'][i, j] + deviation_enec_fjk['wj_fjk_fij_5'][i, j]) == 0
              or (deviation_enec_fjk['wi_fjk_fik_6'][i, j] + deviation_enec_fjk['wk_fjk_fik_6'][i, j]) == 0):
            ENEC_fjk[i, j] = 0
        else: ENEC_fjk[i, j] = 1
print(ENEC_fjk)

situation_initiale = 'fjk'
graph_matrice_e_pythagore(n_b_cij, n_p, titre_b, situation_initiale,deviation_enec_fjk)
graph_matrice_pythagore(n_b_cij, n_p, ENEC_fjk, "ENEC_fjk")


##################################### Déviation de coalition - multilatérale ##############################################

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

ENEC_m_F = np.zeros((n_m_cij, n_p))
deviation_enec_m_F = deviation_coal_matrice_pythagore("mF",symbole_dev, w_F, vect_dev_m_F, vect_cre_m_F, vect_pays, 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
for i in range(n_m_cij):
    for j in range(n_p):
        if EN_m_F[i, j] == 2:
            ENEC_m_F[i, j] = 2
        elif EN_m_F[i, j] == 0:
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
graph_matrice_e_pythagore(n_m_cij, n_p, titre_m, situation_initiale,deviation_enec_m_F)
graph_matrice_pythagore(n_m_cij, n_p, ENEC_m_F, "ENEC_m_F")


# Déviation coalition - statu quo
vect_dev_m_quo = [w_mij, w_mik, w_mjk, w_F]
symbole_dev = ['mij', 'mik', 'mjk', 'mF']
vect_pays = [['i','j'], ['i','k'], ['j','k'], ['i', 'j', 'k']]
vect_cre_m_quo = [[[dev_inverse_mij['wi_quo_mij']],[dev_inverse_mij['wj_quo_mij']]],
                [[dev_inverse_mik['wi_quo_mik']],[dev_inverse_mik['wk_quo_mik']]],
                [[dev_inverse_mjk['wj_quo_mjk']],[dev_inverse_mjk['wk_quo_mjk']]],
                [[inverse_ENEC_m_F], [inverse_ENEC_m_F],[inverse_ENEC_m_F]]
                ]

ENEC_m_quo = np.zeros((n_m_cij, n_p))
deviation_enec_m_quo = deviation_coal_matrice_pythagore('quo', symbole_dev, w_quo, vect_dev_m_quo, vect_cre_m_quo, vect_pays, 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
for i in range(n_m_cij):
    for j in range(n_p):
        if deviation_enec_m_quo['wi_quo_mij_1'][i, j] == 2:
            ENEC_m_quo =2
        elif ((deviation_enec_m_quo['wi_quo_mij_1'][i, j] + deviation_enec_m_quo['wj_quo_mij_1'][i, j]) == 0
              or (deviation_enec_m_quo['wi_quo_mik_2'][i, j] + deviation_enec_m_quo['wk_quo_mik_2'][i, j]) == 0
              or (deviation_enec_m_quo['wj_quo_mjk_3'][i, j] + deviation_enec_m_quo['wk_quo_mjk_3'][i, j]) == 0
              or (deviation_enec_m_quo['wj_quo_mF_4'][i, j] + deviation_enec_m_quo['wk_quo_mF_4'][i, j] + deviation_enec_m_quo['wj_quo_mF_4'][i, j]) == 0):
            ENEC_m_quo[i, j] = 0
        else: ENEC_m_quo[i, j] = 1
print(ENEC_m_quo)

situation_initiale = 'statu quo'
graph_matrice_e_pythagore(n_m_cij, n_p, titre_m, situation_initiale,deviation_enec_m_quo)
graph_matrice_pythagore(n_m_cij, n_p, ENEC_m_quo, "ENEC_m_quo")



# Déviation coalition - mij
dev_manquant_mij, dev_manquant_inverse_mij = deviation_uni_matrice_pythagore("mij",['mik', 'mjk'], w_mij,[w_mik, w_mjk],  ['i','j','k'], 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
dev_manquant_mik, dev_manquant_inverse_mik = deviation_uni_matrice_pythagore("mik",['mij', 'mjk'], w_mik,[w_mij, w_mjk],  ['i','j','k'], 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)

vect_dev_mij = [w_F, w_mik, w_mjk]
symbole_dev = ['mF', 'mik', 'mjk']
vect_pays = [['i', 'j', 'k'], ['i','k'], ['j','k']]
vect_cre_mij = [[[inverse_ENEC_m_F],[inverse_ENEC_m_F], [inverse_ENEC_m_F]],
                [[dev_manquant_mij['wi_mij_mik'], dev_inverse_mik['wi_quo_mik']], [dev_inverse_mik['wk_quo_mik']]],
                [[dev_manquant_mij['wj_mij_mjk'], dev_inverse_mjk['wj_quo_mjk']], [dev_inverse_mjk['wk_quo_mjk']]],
                ]

ENEC_mij = np.zeros((n_m_cij, n_p))
deviation_enec_mij = deviation_coal_matrice_pythagore('mij', symbole_dev, w_mij, vect_dev_mij, vect_cre_mij, vect_pays, 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
for i in range(n_m_cij):
    for j in range(n_p):
        if EN_mij[i, j] == 2:
            ENEC_mij[i, j] = 2
        elif EN_mij[i, j] == 0:
            ENEC_mij[i, j] = 0
        elif ((deviation_enec_mij['wi_mij_mF_1'][i, j] + deviation_enec_mij['wj_mij_mF_1'][i, j] + deviation_enec_mij['wk_mij_mF_1'][i, j]) == 0
              or (deviation_enec_mij['wi_mij_mik_2'][i, j] + deviation_enec_mij['wk_mij_mik_2'][i, j]) == 0
              or (deviation_enec_mij['wj_mij_mjk_3'][i, j] + deviation_enec_mij['wk_mij_mjk_3'][i, j]) == 0):
            ENEC_mij[i, j] = 0
        else: ENEC_mij[i, j] = 1
print(ENEC_mij)

situation_initiale = 'mij'
graph_matrice_e_pythagore(n_m_cij, n_p, titre_m, situation_initiale,deviation_enec_mij)
graph_matrice_pythagore(n_m_cij, n_p, ENEC_mij, "ENEC_mij")


# Déviation coalition - mik
vect_dev_mik = [w_F, w_mij, w_mjk]
symbole_dev = ['mF', 'mij', 'mjk']
vect_pays = [['i', 'j', 'k'], ['i','j'], ['j','k']]
vect_cre_mik = [[[inverse_ENEC_m_F],[inverse_ENEC_m_F], [inverse_ENEC_m_F]],
                [[dev_manquant_mik['wi_mik_mij'], dev_inverse_mij['wi_quo_mij']], [dev_inverse_mij['wj_quo_mij']]],
                [[dev_inverse_mjk['wj_quo_mjk']], [dev_manquant_mik['wk_mik_mjk'], dev_inverse_mjk['wk_quo_mjk']]],
                ]

ENEC_mik = np.zeros((n_m_cij, n_p))
deviation_enec_mik = deviation_coal_matrice_pythagore('mik', symbole_dev, w_mik, vect_dev_mik, vect_cre_mik, vect_pays, 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
for i in range(n_m_cij):
    for j in range(n_p):
        if EN_mik[i, j] == 2:
            ENEC_mik[i, j] = 2
        elif EN_mik[i, j] == 0:
            ENEC_mik[i, j] = 0
        elif ((deviation_enec_mik['wi_mik_mF_1'][i, j] + deviation_enec_mik['wj_mik_mF_1'][i, j] + deviation_enec_mik['wk_mik_mF_1'][i, j]) == 0
              or (deviation_enec_mik['wi_mik_mij_2'][i, j] + deviation_enec_mik['wj_mik_mij_2'][i, j]) == 0
              or (deviation_enec_mik['wj_mik_mjk_3'][i, j] + deviation_enec_mik['wk_mik_mjk_3'][i, j]) == 0):
            ENEC_mik[i, j] = 0
        else: ENEC_mik[i, j] = 1
print(ENEC_mik)

situation_initiale = 'mik'
graph_matrice_e_pythagore(n_m_cij, n_p, titre_m, situation_initiale,deviation_enec_mik)
graph_matrice_pythagore(n_m_cij, n_p, ENEC_mik, "ENEC_mik")



# Déviation coalition - mjk
vect_dev_mjk = [w_F, w_mij, w_mik]
symbole_dev = ['mF', 'mij', 'mik']
vect_pays = [['i', 'j', 'k'], ['i','j'], ['i','k']]
vect_cre_mjk = [[[inverse_ENEC_m_F],[inverse_ENEC_m_F], [inverse_ENEC_m_F]],
                [[dev_inverse_mij['wi_quo_mij']], [dev_manquant_inverse_mij['wj_mjk_mij'], dev_inverse_mij['wj_quo_mij']]],
                [[dev_inverse_mik['wi_quo_mik']], [dev_manquant_inverse_mik['wk_mjk_mik'], dev_inverse_mik['wk_quo_mik']]],
                ]

ENEC_mjk = np.zeros((n_m_cij, n_p))
deviation_enec_mjk = deviation_coal_matrice_pythagore('mjk', symbole_dev, w_mjk, vect_dev_mjk, vect_cre_mjk, vect_pays, 'multilatérale', n_m_cij, n_p, c_values_m_cij, p_values)
for i in range(n_m_cij):
    for j in range(n_p):
        if EN_mjk[i, j] == 2:
            ENEC_mjk[i, j] = 2
        elif EN_mjk[i, j] == 0:
            ENEC_mjk[i, j] = 0
        elif ((deviation_enec_mjk['wi_mjk_mF_1'][i, j] + deviation_enec_mjk['wj_mjk_mF_1'][i, j] + deviation_enec_mjk['wk_mjk_mF_1'][i, j]) == 0
              or (deviation_enec_mjk['wi_mjk_mij_2'][i, j] + deviation_enec_mjk['wj_mjk_mij_2'][i, j]) == 0
              or (deviation_enec_mjk['wi_mjk_mik_3'][i, j] + deviation_enec_mjk['wk_mjk_mik_3'][i, j]) == 0):
            ENEC_mjk[i, j] = 0
        else: ENEC_mjk[i, j] = 1
print(ENEC_mjk)

situation_initiale = 'mjk'
graph_matrice_e_pythagore(n_m_cij, n_p, titre_m, situation_initiale,deviation_enec_mjk)
graph_matrice_pythagore(n_m_cij, n_p, ENEC_mjk, "ENEC_mjk")

