import numpy as np
from sympy import *
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sympy import solve, Rational,simplify
from sympy.conftest import sp
from sympy.core.relational import Relational
ei, ej, ek, cij, cik, cji, cjk, cki, ckj, tij, tji, tik, tki, tjk, tkj = symbols('ei ej ek cij cik cji cjk cki ckj '
                                                                                 'tij tji tik tki tjk tkj')
e, c, a, p = symbols('e c a p')
ci, cj = symbols('ci cj')

##############################################################################################################
##############################################################################################################
##############################################################################################################
########################################## Fonctions générales ###############################################

def prix_equilibre_1(a, ei, ej, ek, ti_quo, tj_quo, tk_quo, ti_fij, tj_fij, ti_fik, tk_fik, tj_fjk, tk_fjk,ti_mij,
                   tj_mij, ti_mik, tk_mik, tj_mjk, tk_mjk, tj_ih, tk_ih, ti_jh, tk_jh, ti_kh, tj_kh, cij, cji, cik,
                     cki, cjk, ckj):
    # fjk
    pi_I_fjk = (Rational(1, 3) * (3 * a - (ej + ek) + (ti_quo + ti_quo) + (cij + cik)))
    pi_J_fjk = (Rational(1, 3) * (3 * a - (ei + ek) - 2*tj_fjk + 0 - 2*cji + cjk))
    pi_K_fjk = (Rational(1, 3) * (3 * a - (ei + ej) - 2*tk_fjk + 0 - 2*cki + ckj))

    pj_J_fjk = (Rational(1, 3) * (3 * a - (ei + ek) + (tj_fjk + 0) + (cji + cjk)))
    pj_I_fjk = (Rational(1, 3) * (3 * a - (ej + ek) - 2 * ti_quo + ti_quo - 2 * cij + cik))
    pj_K_fjk = (Rational(1, 3) * (3 * a - (ei + ej) - 0 + tk_fjk - 2 * ckj + cki))

    pk_K_fjk = (Rational(1, 3) * (3 * a - (ej + ei) + (0 + tk_fjk) + (ckj + cki)))
    pk_I_fjk = (Rational(1, 3) * (3 * a - (ej + ek) - 2 * ti_quo + ti_quo - 2 * cik + cij))
    pk_J_fjk = (Rational(1, 3) * (3 * a - (ei + ek) - 0 + tj_fjk - 2 * cjk + cji))

    # mjk
    pi_I_mjk = (Rational(1, 3) * (3 * a - (ej + ek) + (ti_quo + ti_quo) + (cij + cik)))
    pi_J_mjk = (Rational(1, 3) * (3 * a - (ei + ek) - 2 * tj_mjk + tj_mjk - 2 * cji + cjk))
    pi_K_mjk = (Rational(1, 3) * (3 * a - (ei + ej) - 2 * tk_mjk + tk_mjk - 2 * cki + ckj))

    pj_J_mjk = (Rational(1, 3) * (3 * a - (ei + ek) + (tj_mjk + tj_mjk) + (cji + cjk)))
    pj_I_mjk = (Rational(1, 3) * (3 * a - (ej + ek) - 2 * ti_quo + ti_quo - 2 * cij + cik))
    pj_K_mjk = (Rational(1, 3) * (3 * a - (ei + ej) - 2 * tk_mjk + tk_mjk - 2 * ckj + cki))

    pk_K_mjk = (Rational(1, 3) * (3 * a - (ej + ei) + (tk_mjk + tk_mjk) + (ckj + cki)))
    pk_I_mjk = (Rational(1, 3) * (3 * a - (ej + ek) - 2 * ti_quo + ti_quo - 2 * cij + cik))
    pk_J_mjk = (Rational(1, 3) * (3 * a - (ei + ek) - 2 * tj_mjk + tj_mjk - 2 * ckj + cki))



    P_vecteur_1 = np.array([pi_I_fjk, pi_J_fjk, pi_K_fjk, pj_J_fjk, pj_I_fjk, pj_K_fjk, pk_K_fjk, pk_I_fjk, pk_J_fjk,
                          pi_I_mjk, pi_J_mjk, pi_K_mjk, pj_J_mjk, pj_I_mjk, pj_K_mjk, pk_K_mjk, pk_I_mjk, pk_J_mjk])
    P_dict_1 = {
        "pi_I_fjk": pi_I_fjk, "pi_J_fjk": pi_J_fjk, "pi_K_fjk": pi_K_fjk,
        "pj_J_fjk": pj_J_fjk, "pj_I_fjk": pj_I_fjk, "pj_K_fjk": pj_K_fjk,
        "pk_K_fjk": pk_K_fjk, "pk_I_fjk": pk_I_fjk, "pk_J_fjk": pk_J_fjk,
        "pi_I_mjk": pi_I_mjk, "pi_J_mjk": pi_J_mjk, "pi_K_mjk": pi_K_mjk,
        "pj_J_mjk": pj_J_mjk, "pj_I_mjk": pj_I_mjk, "pj_K_mjk": pj_K_mjk,
        "pk_K_mjk": pk_K_mjk, "pk_I_mjk": pk_I_mjk, "pk_J_mjk": pk_J_mjk,
    }
    for key, i in zip(P_dict_1,range(len(P_vecteur_1))):  # Utilisation de len(T) pour éviter les erreurs d'index
        print(f"[{key}] = {simplify(P_vecteur_1[i])}\n")
    return P_vecteur_1, P_dict_1


def prix_equilibre(a, ei, ej, ek, ti_quo, tj_quo, tk_quo, ti_fij, tj_fij, ti_fik, tk_fik, tj_fjk, tk_fjk,ti_mij,
                   tj_mij, ti_mik, tk_mik, tj_mjk, tk_mjk, tj_ih, tk_ih, ti_jh, tk_jh, ti_kh, tj_kh, cij, cji, cik,
                   cki, cjk, ckj):
    #quo
    pi_i_quo = (Rational(1, 3)*(3*a - (ej+ek) + (ti_quo+ti_quo) + (cij+cik)))

    pj_j_quo = (Rational(1, 3)*(3*a - (ei+ek) + (tj_quo+tj_quo) + (cji+cjk)))

    pk_k_quo = (Rational(1, 3)*(3*a - (ej+ei) + (tk_quo+tk_quo) + (ckj+cki)))

    # fij
    pi_i_fij = (Rational(1, 3) * (3 * a - (ej + ek) + (0 + ti_fij) + (cij + cik)))

    pj_j_fij = (Rational(1, 3) * (3 * a - (ei + ek) + (0 + tj_fij) + (cji + cjk)))

    pk_k_fij = (Rational(1, 3) * (3 * a - (ej + ei) + (tk_quo + tk_quo) + (ckj + cki)))

    # fik
    pi_i_fik = (Rational(1, 3) * (3 * a - (ej + ek) + (ti_fik + 0) + (cij + cik)))

    pj_j_fik = (Rational(1, 3) * (3 * a - (ei + ek) + (tj_quo + tj_quo) + (cji + cjk)))

    pk_k_fik = (Rational(1, 3) * (3 * a - (ej + ei) + (tk_fik + 0) + (ckj + cki)))

    # fjk
    pi_i_fjk = (Rational(1, 3) * (3 * a - (ej + ek) + (ti_quo + ti_quo) + (cij + cik)))

    pj_j_fjk = (Rational(1, 3) * (3 * a - (ei + ek) + (tj_fjk + 0) + (cji + cjk)))

    pk_k_fjk = (Rational(1, 3) * (3 * a - (ej + ei) + (0 + tk_fjk) + (ckj + cki)))

    # mij
    pi_i_mij = (Rational(1, 3) * (3 * a - (ej + ek) + (ti_mij + ti_mij) + (cij + cik)))

    pj_j_mij = (Rational(1, 3) * (3 * a - (ei + ek) + (tj_mij + tj_mij) + (cji + cjk)))

    pk_k_mij = (Rational(1, 3) * (3 * a - (ej + ei) + (tk_quo + tk_quo) + (ckj + cki)))

    # mik
    pi_i_mik = (Rational(1, 3) * (3 * a - (ej + ek) + (ti_mik + ti_mik) + (cij + cik)))

    pj_j_mik = (Rational(1, 3) * (3 * a - (ei + ek) + (tj_quo + tj_quo) + (cji + cjk)))

    pk_k_mik = (Rational(1, 3) * (3 * a - (ej + ei) + (tk_mik + tk_mik) + (ckj + cki)))

    # mjk
    pi_i_mjk = (Rational(1, 3) * (3 * a - (ej + ek) + (ti_quo + ti_quo) + (cij + cik)))

    pj_j_mjk = (Rational(1, 3) * (3 * a - (ei + ek) + (tj_mjk + tj_mjk) + (cji + cjk)))

    pk_k_mjk = (Rational(1, 3) * (3 * a - (ej + ei) + (tk_mjk + tk_mjk) + (ckj + cki)))

    # ih
    pi_i_ih = (Rational(1, 3) * (3 * a - (ej + ek) + (0 + 0) + (cij + cik)))

    pj_j_ih = (Rational(1, 3) * (3 * a - (ei + ek) + (0 + tj_ih) + (cji + cjk)))

    pk_k_ih = (Rational(1, 3) * (3 * a - (ej + ei) + (tk_ih + 0) + (ckj + cki)))

    # jh
    pi_i_jh = (Rational(1, 3) * (3 * a - (ej + ek) + (0 + ti_jh) + (cij + cik)))

    pj_j_jh = (Rational(1, 3) * (3 * a - (ei + ek) + (0 + 0) + (cji + cjk)))

    pk_k_jh = (Rational(1, 3) * (3 * a - (ej + ei) + (0 + tk_jh) + (ckj + cki)))

    # kh
    pi_i_kh = (Rational(1, 3) * (3 * a - (ej + ek) + (ti_kh + 0) + (cij + cik)))

    pj_j_kh = (Rational(1, 3) * (3 * a - (ei + ek) + (tj_kh + 0) + (cji + cjk)))

    pk_k_kh = (Rational(1, 3) * (3 * a - (ej + ei) + (0 + 0) + (ckj + cki)))

    # F
    pi_i_F = (Rational(1, 3) * (3 * a - (ej + ek) + (0 + 0) + (cij + cik)))

    pj_j_F = (Rational(1, 3) * (3 * a - (ei + ek) + (0 + 0) + (cji + cjk)))

    pk_k_F = (Rational(1, 3) * (3 * a - (ej + ei) + (0 + 0) + (ckj + cki)))

    P_vecteur = np.array([pi_i_quo, pj_j_quo, pk_k_quo, pi_i_fij, pj_j_fij, pk_k_fij, pi_i_fik, pj_j_fik, pk_k_fik,
                          pi_i_fjk, pj_j_fjk, pk_k_fjk
                  , pi_i_mij, pj_j_mij, pk_k_mij, pi_i_mik, pj_j_mik, pk_k_mik, pi_i_mjk, pj_j_mjk, pk_k_mjk, pi_i_ih,
                          pj_j_ih, pk_k_ih
                  , pi_i_jh, pj_j_jh, pk_k_jh, pi_i_kh, pj_j_kh, pk_k_kh, pi_i_F, pj_j_F, pk_k_F])
    P_dict = {
        "pi_i_quo": pi_i_quo, "pj_j_quo": pj_j_quo, "pk_k_quo": pk_k_quo,
        "pi_i_fij": pi_i_fij, "pj_j_fij": pj_j_fij, "pk_k_fij": pk_k_fij,
        "pi_i_fik": pi_i_fik, "pj_j_fik": pj_j_fik, "pk_k_fik": pk_k_fik,
        "pi_i_fjk": pi_i_fjk, "pj_j_fjk": pj_j_fjk, "pk_k_fjk": pk_k_fjk,
        "pi_i_mij": pi_i_mij, "pj_j_mij": pj_j_mij, "pk_k_mij": pk_k_mij,
        "pi_i_mik": pi_i_mik, "pj_j_mik": pj_j_mik, "pk_k_mik": pk_k_mik,
        "pi_i_mjk": pi_i_mjk, "pj_j_mjk": pj_j_mjk, "pk_k_mjk": pk_k_mjk,
        "pi_i_ih": pi_i_ih, "pj_j_ih": pj_j_ih, "pk_k_ih": pk_k_ih,
        "pi_i_jh": pi_i_jh, "pj_j_jh": pj_j_jh, "pk_k_jh": pk_k_jh,
        "pi_i_kh": pi_i_kh, "pj_j_kh": pj_j_kh, "pk_k_kh": pk_k_kh,
        "pi_i_F": pi_i_F, "pj_j_F": pj_j_F, "pk_k_F": pk_k_F
    }

    return P_vecteur, P_dict


def export_i(a, ei, ej, ek, pi_i, pj_j, pk_k, tij, tik, tji, tjk, tki, tkj, cij, cji, cik, cki, cjk, ckj):
    XI_j = ej - a + pi_i - tij - cij

    XI_k = ek - a + pi_i - tik - cik

    XJ_i = ei - a + pj_j - tji - cji

    XJ_k = ek - a + pj_j - tjk - cjk

    XK_i = ei - a + pk_k - tki - cki

    XK_j = ej - a + pk_k - tkj - ckj

    X = [XI_j, XI_k, XJ_i, XJ_k, XK_i, XK_j]

    return X
def export_final(x_quo, x_fij, x_fik, x_fjk, x_mij, x_mik, x_mjk, x_ih, x_jh, x_kh, x_F):
    X_b = []
    X_b.extend(x_quo)
    X_b.extend(x_fij)
    X_b.extend(x_fik)
    X_b.extend(x_fjk)
    X_b.extend(x_ih)
    X_b.extend(x_jh)
    X_b.extend(x_kh)
    X_b.extend(x_F)

    X_m = []
    X_m.extend(x_mij)
    X_m.extend(x_mik)
    X_m.extend(x_mjk)
    return X_b, X_m



def tarifs_optimaux_bilateraux(ei, ej, ek, cij, cji, cik, cki, cjk, ckj):
    ti_quo = ((ej + ek) - (cij + cik)) / 8  # tij=tik=ti_quo
    tj_quo = ((ei + ek) - (cji + cjk)) / 8  # tji=tjk=tj_quo
    tk_quo = ((ei + ej) - (cki + ckj)) / 8  # tki=tkj=tk_quo

    ti_fij = ((5 * ek - 4 * ej) + (4 * cij - 5 * cik)) / 11  # tij=tji=0 , tik=ti_fij, tki=tkj=tk_quo
    tj_fij = ((5 * ek - 4 * ei) + (4 * cji - 5 * cjk)) / 11  # tjk=tj_fij

    ti_fik = ((5 * ej - 4 * ek) + (4 * cik - 5 * cij)) / 11  # tik=tki=0 , tij=ti_fik, tji=tjk=tj_quo
    tk_fik = ((5 * ej - 4 * ei) + (4 * cki - 5 * ckj)) / 11  # tkj=tk_fik

    tj_fjk = ((5 * ei - 4 * ek) + (4 * cjk - 5 * cji)) / 11  # tjk=tkj=0 , tji=tj_fjk , tij=ti_quo
    tk_fjk = ((5 * ei - 4 * ej) + (4 * ckj - 5 * cki)) / 11  # tki=tk_fjk , tik=ti_quo

    tj_ih = ((5 * ek - 4 * ei) + (4 * cji - 5 * cjk)) / 11  # tij=tik=0 , tji=tki=0 , tjk = tj_ih
    tk_ih = ((5 * ej - 4 * ei) + (4 * cki - 5 * ckj)) / 11  # tkj = tk_ih

    ti_jh = ((5 * ek - 4 * ej) + (4 * cij - 5 * cik)) / 11  # tij=tkj=0 , tji=tjk=0, tik = ti_jh
    tk_jh = ((5 * ei - 4 * ej) + (4 * ckj - 5 * cki)) / 11  # tki = tk_jh

    ti_kh = ((5 * ej - 4 * ek) + (4 * cik - 5 * cij)) / 11  # tki=tkj=0 , tik=tjk=0, tij = ti_kh
    tj_kh = ((5 * ei - 4 * ek) + (4 * cjk - 5 * cji)) / 11  # tji = tj_kh

    T_bilateraux_vecteur = [ti_quo, tj_quo, tk_quo, ti_fij, tj_fij, ti_fik, tk_fik, tj_fjk, tk_fjk, tj_ih, tk_ih,
                            ti_jh, tk_jh, ti_kh, tj_kh]

    T_bilateraux_dict = {
        "ti_quo": ti_quo,
        "tj_quo": tj_quo,
        "tk_quo": tk_quo,
        "ti_fij": ti_fij,
        "tj_fij": tj_fij,
        "ti_fik": ti_fik,
        "tk_fik": tk_fik,
        "tj_fjk": tj_fjk,
        "tk_fjk": tk_fjk,
        "tj_ih": tj_ih,
        "tk_ih": tk_ih,
        "ti_jh": ti_jh,
        "tk_jh": tk_jh,
        "ti_kh": ti_kh,
        "tj_kh": tj_kh
    }
    for key, i in zip(T_bilateraux_dict, range(len(T_bilateraux_vecteur))):
        print(f"T_bilateraux[{key}] = {simplify(T_bilateraux_vecteur[i])}\n")
    print("")
    return T_bilateraux_vecteur, T_bilateraux_dict

def tarifs_optimaux_multilateraux(ei, ej, ek, cij, cji, cik, cki, cjk, ckj):
    ti_quo = ((ej + ek) - (cij + cik)) / 8  # tij=tik=ti_quo
    tj_quo = ((ei + ek) - (cji + cjk)) / 8  # tji=tjk=tj_quo
    tk_quo = ((ei + ej) - (cki + ckj)) / 8  # tki=tkj=tk_quo

    ti_mij = ((2 * ek - ej) + (cij - 2 * cik)) / 7  # tij=tik=ti_mij , tki=tkj=tk_quo
    tj_mij = ((2 * ek - ei) + (cji - 2 * cjk)) / 7  # tji=tjk=tj_mij

    ti_mik = ((2 * ej - ek) + (cik - 2 * cij)) / 7  # tij=tik=ti_mij , tki=tkj=tk_quo
    tk_mik = ((2 * ej - ei) + (cki - 2 * ckj)) / 7  # tji=tjk=tj_mij

    tj_mjk = ((2 * ei - ek) + (cjk - 2 * cji)) / 7  # tji=tjk=tj_mjk , tij=tik-ti_quo
    tk_mjk = ((2 * ei - ej) + (ckj - 2 * cki)) / 7  # tki=tkj=tk_mjk

    T_multilateraux_vecteur = [ti_quo, tj_quo, tk_quo, ti_mij, tj_mij,ti_mik, tk_mik, tj_mjk, tk_mjk]
    T_multilateraux_dict = {
        "ti_quo": ti_quo,
        "tj_quo": tj_quo,
        "tk_quo": tk_quo,
        "ti_mij": ti_mij,
        "tj_mij": tj_mij,
        "ti_mik": ti_mik,
        "tk_mik": tk_mik,
        "tj_mjk": tj_mjk,
        "tk_mjk": tk_mjk,
    }
    for key, i in zip(T_multilateraux_dict, range(len(T_multilateraux_vecteur))):
        print(f"T_multilateraux[{key}] = {simplify(T_multilateraux_vecteur[i])}\n")
    print("")
    return T_multilateraux_vecteur, T_multilateraux_dict



def w_optimaux_sep(a, ei, ej, ek, cij, cji, cik, cki, cjk, ckj, tij, tji, tik, tki, tjk, tkj, regime):
    sci = (Rational(1, 18) * (
                ((ej + ek) - (tij + tik) - (cij + cik)) ** 2 + ((ei + ek) + (2 * tji - tjk) + (2 * cji - cjk)) ** 2
                + ((ei + ej) + (2 * tki - tkj) + (2 * cki - ckj)) ** 2))

    spi = (Rational(1, 3) * (ei) * (6 * a - (ei + ek) + (tji + tjk) + (cji + cjk) - 3 * tji - 3 * cji - (ei + ej)
                                    + (tki + tkj) +
                      (cki + ckj) - 3 * tki - 3 * cki))
    gri = (Rational(1, 3) * tij * ((2 * ej - ek) + (tik - 2 * tij) + (cik - 2 * cij)) + Rational(1, 3) *
           tik * (
                      (2 * ek - ej) + (tij - 2 * tik) + (cij - 2 * cik)))

    scj = (Rational(1, 18) * (
                ((ei + ek) - (tji + tjk) - (cji + cjk)) ** 2 + ((ej + ek) + (2 * tij - tik) + (2 * cij - cik)) ** 2
                + ((ej + ei) + (2 * tkj - tki) + (2 * ckj - cki)) ** 2))
    spj = (Rational(1, 3) * (ej) * (
                      6 * a - (ej + ek) + (tij + tik) + (cij + cik) - 3 * tij - 3 * cij - (ej + ei) + (tkj + tki) +
                      (ckj + cki) - 3 * tkj - 3 * ckj))
    grj = (Rational(1, 3) * tji * ((2 * ei - ek) + (tjk - 2 * tji) + (cjk - 2 * cji)) + Rational(1, 3) *
           tjk * (
                      (2 * ek - ei) + (tji - 2 * tjk) + (cji - 2 * cjk)))

    sck = (Rational(1, 18) * (
                ((ej + ei) - (tkj + tki) - (ckj + cki)) ** 2 + ((ek + ei) + (2 * tjk - tji) + (2 * cjk - cji)) ** 2
                + ((ek + ej) + (2 * tik - tij) + (2 * cik - cij)) ** 2))
    spk = ( Rational(1, 3) * (ek) * (
                      6 * a - (ek + ei) + (tjk + tji) + (cjk + cji) - 3 * tjk - 3 * cjk - (ek + ej) + (tik + tij) +
                      (cik + cij) - 3 * tik - 3 * cik))
    grk = (Rational(1, 3) * tkj * ((2 * ej - ei) + (tki - 2 * tkj) + (cki - 2 * ckj)) + Rational(1, 3) *
           tki * (
                      (2 * ei - ej) + (tkj - 2 * tki) + (ckj - 2 * cki)))
    sc_vect = np.array([sci, scj, sck])
    sp_vect = np.array([spi, spj, spk])
    gr_vect = np.array([gri, grj, grk])
    sc_dict = {
        f"sci_{regime}": sci,
        f"scj_{regime}": scj,
        f"sck_{regime}": sck
    }
    sp_dict = {
        f"spi_{regime}": spi,
        f"spj_{regime}": spj,
        f"spk_{regime}": spk
    }
    gr_dict = {
        f"gri_{regime}": gri,
        f"grj_{regime}": grj,
        f"grk_{regime}": grk
    }
    for key, i in zip(sc_dict, range(len(sc_vect))):
        print(f"{key} = {(sc_vect[i])}\n")
    for key, i in zip(sp_dict, range(len(sp_vect))):
        print(f"{key} = {(sp_vect[i])}\n")
    for key, i in zip(gr_dict, range(len(gr_vect))):
        print(f"{key} = {(gr_vect[i])}\n")

    return


def w_optimaux(a, ei, ej, ek, cij, cji, cik, cki, cjk, ckj, tij, tji, tik, tki, tjk, tkj, regime):
    wi = (Rational(1, 18) * (
                ((ej + ek) - (tij + tik) - (cij + cik)) ** 2 + ((ei + ek) + (2 * tji - tjk) + (2 * cji - cjk)) ** 2
                + ((ei + ej) + (2 * tki - tkj) + (2 * cki - ckj)) ** 2)
          + Rational(1, 3) * (ei) * (
                      6 * a - (ei + ek) + (tji + tjk) + (cji + cjk) - 3 * tji - 3 * cji - (ei + ej) + (tki + tkj) +
                      (cki + ckj) - 3 * tki - 3 * cki)
          + Rational(1, 3) * tij * ((2 * ej - ek) + (tik - 2 * tij) + (cik - 2 * cij)) + Rational(1, 3) *
          tik * (
                      (2 * ek - ej) + (tij - 2 * tik) + (cij - 2 * cik)))

    wj = (Rational(1, 18) * (
                ((ei + ek) - (tji + tjk) - (cji + cjk)) ** 2 + ((ej + ek) + (2 * tij - tik) + (2 * cij - cik)) ** 2
                + ((ej + ei) + (2 * tkj - tki) + (2 * ckj - cki)) ** 2)
          + Rational(1, 3) * (ej) * (
                      6 * a - (ej + ek) + (tij + tik) + (cij + cik) - 3 * tij - 3 * cij - (ej + ei) + (tkj + tki) +
                      (ckj + cki) - 3 * tkj - 3 * ckj)
          + Rational(1, 3) * tji * ((2 * ei - ek) + (tjk - 2 * tji) + (cjk - 2 * cji)) + Rational(1, 3) *
          tjk * (
                      (2 * ek - ei) + (tji - 2 * tjk) + (cji - 2 * cjk)))

    wk = (Rational(1, 18) * (
                ((ej + ei) - (tkj + tki) - (ckj + cki)) ** 2 + ((ek + ei) + (2 * tjk - tji) + (2 * cjk - cji)) ** 2
                + ((ek + ej) + (2 * tik - tij) + (2 * cik - cij)) ** 2)
          + Rational(1, 3) * (ek) * (
                      6 * a - (ek + ei) + (tjk + tji) + (cjk + cji) - 3 * tjk - 3 * cjk - (ek + ej) + (tik + tij) +
                      (cik + cij) - 3 * tik - 3 * cik)
          + Rational(1, 3) * tkj * ((2 * ej - ei) + (tki - 2 * tkj) + (cki - 2 * ckj)) + Rational(1, 3) *
          tki * (
                      (2 * ei - ej) + (tkj - 2 * tki) + (ckj - 2 * cki)))
    W_vect = np.array([wi, wj, wk])
    W_dict = {
        f"wi_{regime}": wi,
        f"wj_{regime}": wj,
        f"wk_{regime}": wk
    }
    for key, i in zip(W_dict, range(len(W_vect))):
        print(f"{key} = {simplify(W_vect[i])}\n")

    return W_vect
def w_optimaux_virgule(a, ei, ej, ek, cij, cji, cik, cki, cjk, ckj, tij, tji, tik, tki, tjk, tkj, regime):
    wi = ((1/18) * (((ej+ek) - (tij+tik) - (cij+cik))**2 + ((ei+ek) + (2*tji-tjk) + (2*cji-cjk))**2
                     + ((ei+ej) + (2*tki-tkj) + (2*cki-ckj))**2)
          + (ei/3) * (6*a - (ei+ek) + (tji+tjk) + (cji+cjk) - 3*tji - 3*cji - (ei+ej) + (tki+tkj) +
                        (cki+ckj) - 3*tki - 3*cki)
          + (1/3)*tij * ((2*ej-ek) + (tik-2*tij) + (cik-2*cij)) + (1/3)*tik * ((2*ek-ej) + (tij-2*tik) + (cij-2*cik)))

    wj = ((1/18) * (((ei+ek) - (tji+tjk) - (cji+cjk))**2 + ((ej+ek) + (2*tij-tik) + (2*cij-cik))**2
                     + ((ej+ei) + (2*tkj-tki) + (2*ckj-cki))**2)
          + (ej/3) * (6*a - (ej+ek) + (tij+tik) + (cij+cik) - 3*tij - 3*cij - (ej+ei) + (tkj+tki) +
                        (ckj+cki) - 3*tkj - 3*ckj)
          + (1/3)*tji * ((2*ei-ek) + (tjk-2*tji) + (cjk-2*cji)) + (1/3)*tjk * ((2*ek-ei) + (tji-2*tjk) + (cji-2*cjk)))

    wk = ((1/18) * (((ej+ei) - (tkj+tki) - (ckj+cki))**2 + ((ek+ei) + (2*tjk-tji) + (2*cjk-cji))**2
                     + ((ek+ej) + (2*tik-tij) + (2*cik-cij))**2)
          + (ek/3) * (6*a - (ek+ei) + (tjk+tji) + (cjk+cji) - 3*tjk - 3*cjk - (ek+ej) + (tik+tij) +
                        (cik+cij) - 3*tik - 3*cik)
          + (1/3)*tkj * ((2*ej-ei) + (tki-2*tkj) + (cki-2*ckj)) + (1/3)*tki * ((2*ei-ej) + (tkj-2*tki) + (ckj-2*cki)))
    W_vect = np.array([wi, wj, wk])
    W_dict = {
        f"wi_{regime}": wi,
        f"wj_{regime}": wj,
        f"wk_{regime}": wk
    }
    for key, i in zip(W_dict, range(len(W_vect))):
        print(f"{key} = {simplify((W_vect[i]))}\n")

    return W_vect


##############################################################################################################
##############################################################################################################
############################################# Situation 1 ####################################################
def obtenir_bornes_restrictives(s_t_o_b_i, s_t_o_b_i_dict, s_t_o_m_i, s_t_o_m_i_dict, X_b, X_m, variable):
    def normaliser_inegalites(inegalites):
        bornes = []
        for ineq in inegalites:
            if isinstance(ineq, Relational):  # Vérifier si c'est une inégalité
                if ineq.lhs == variable:  # Exemple: c <= 2 ou c >= -1
                    bornes.append(ineq) # append, car un seul élément
                elif ineq.rhs == variable:  # Exemple: -2/3 <= c  devient  c >= -2/3
                    bornes.append(ineq.reversed)
        return bornes
    bornes_bilaterales = []
    bornes_multilaterales = []
    bornes_b_max = []
    bornes_b_min = []
    bornes_m_max = []
    bornes_m_min = []
    bornes_restrictives_bilaterale = []
    bornes_restrictives_multilaterale = []
    for tarif in range(len(s_t_o_b_i)):
        result = solve(s_t_o_b_i[tarif]>=0, variable)
        bornes_bilaterales.extend(result.args) # .extend() , car solve() renvoit 2 conditions en même temps
    for tarif in range(len(s_t_o_m_i)):
        result = solve(s_t_o_m_i[tarif]>=0, variable)
        bornes_multilaterales.extend(result.args)
    for x in range(len(X_b)):
        result = solve(X_b[x]>=0, variable)
        bornes_bilaterales.extend(result.args)
    for x in range(len(X_m)):
        result = solve(X_m[x]>=0, variable)
        bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["ti_quo"]-s_t_o_b_i_dict["ti_fij"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_fik"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_jh"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_kh"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_fij"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_fjk"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_ih"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_kh"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_fik"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_fjk"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_ih"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_jh"] >= 0, variable)
    bornes_bilaterales.extend(result.args)
    result = solve(s_t_o_m_i_dict["ti_quo"] - s_t_o_m_i_dict["ti_mij"] >= 0, variable)
    bornes_multilaterales.extend(result.args)
    result = solve(s_t_o_m_i_dict["ti_quo"] - s_t_o_m_i_dict["ti_mik"] >= 0, variable)
    bornes_multilaterales.extend(result.args)
    result = solve(s_t_o_m_i_dict["tj_quo"] - s_t_o_m_i_dict["tj_mij"] >= 0, variable)
    bornes_multilaterales.extend(result.args)
    result = solve(s_t_o_m_i_dict["tj_quo"] - s_t_o_m_i_dict["tj_mjk"] >= 0, variable)
    bornes_multilaterales.extend(result.args)
    result = solve(s_t_o_m_i_dict["tk_quo"] - s_t_o_m_i_dict["tk_mik"] >= 0, variable)
    bornes_multilaterales.extend(result.args)
    result = solve(s_t_o_m_i_dict["tk_quo"] - s_t_o_m_i_dict["tk_mjk"] >= 0, variable)
    bornes_multilaterales.extend(result.args)
    bornes_bilaterales = normaliser_inegalites(bornes_bilaterales)
    bornes_multilaterales = normaliser_inegalites(bornes_multilaterales)

    for borne in bornes_bilaterales:
        if isinstance(borne, Relational) and borne.rel_op == '<=':
            valeur = borne.rhs
            bornes_b_max.append(valeur)
    print("")

    for borne in bornes_bilaterales:
        if isinstance(borne, Relational) and borne.rel_op == '>=':
            valeur = borne.rhs
            bornes_b_min.append(valeur)

    bornes_restrictives_bilaterale.append(Max(*bornes_b_min) if Max(*bornes_b_min) > 0 else 0)
    bornes_restrictives_bilaterale.append(Min(*bornes_b_max))

    for borne in bornes_multilaterales:
        if isinstance(borne, Relational) and borne.rel_op == '<=':
            valeur = borne.rhs
            bornes_m_max.append(valeur)
    print("")

    for borne in bornes_multilaterales:
        if isinstance(borne, Relational) and borne.rel_op == '>=':
            valeur = borne.rhs
            bornes_m_min.append(valeur)

    bornes_restrictives_multilaterale.append(Max(*bornes_m_min) if Max(*bornes_m_min) > 0 else 0)
    bornes_restrictives_multilaterale.append(Min(*bornes_m_max))
    print(f'borne restrictive bilatérale = {bornes_restrictives_bilaterale}\n')
    print(f'borne restrictive multilatérale = {bornes_restrictives_multilaterale}\n')
    return bornes_bilaterales, bornes_multilaterales, bornes_restrictives_bilaterale, bornes_restrictives_multilaterale


def axe(structure, vecteur, ax, c_values, position):
    if structure == "bilatérale":
        if position == 0:
            if all(i == 0 for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['red']),
                           marker="s", label=f"{vecteur}")
            elif all(i == 1 for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['#00FF00']),
                           marker="s", label=f"{vecteur}")
            elif all(i in {0, 2} for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['red', 'black']),
                           marker="s", label=f"{vecteur}")
            elif all(i in {1, 2} for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['#00FF00', 'black']),
                           marker="s", label=f"{vecteur}")
            elif all(i in {0, 1, 2} for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['red', '#00FF00', 'black']),
                           marker="s", label=f"{vecteur}")
            elif all(i == 3 for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['white']),
                           marker="s", label=f"{vecteur}")
            else:
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['red', '#00FF00']),
                           marker="s", label=f"{vecteur}")
        else:
            if all(i == 0 for i in vecteur):
                ax.scatter(c_values, position * np.ones_like(c_values), c=vecteur, cmap=mcolors.ListedColormap(['red']),
                           marker="s", label=f"{vecteur}")
            elif all(i == 1 for i in vecteur):
                ax.scatter(c_values, position * np.ones_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['#00FF00']),
                           marker="s", label=f"{vecteur}")
            elif all(i in {0, 2} for i in vecteur):
                ax.scatter(c_values, position * np.ones_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['red', 'black']),
                           marker="s", label=f"{vecteur}")
            elif all(i in {1, 2} for i in vecteur):
                ax.scatter(c_values, position * np.ones_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['#00FF00', 'black']),
                           marker="s", label=f"{vecteur}")
            elif all(i in {0, 1, 2} for i in vecteur):
                ax.scatter(c_values, position * np.ones_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['red', '#00FF00', 'black']),
                           marker="s", label=f"{vecteur}")
            elif all(i == 3 for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['white']),
                           marker="s", label=f"{vecteur}")
            else:
                ax.scatter(c_values, position * np.ones_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['red', '#00FF00']),
                           marker="s", label=f"{vecteur}")
    else:
        if position == 0:
            if all(i == 0 for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['red']),
                           marker="s", label=f"{vecteur}")
            elif all(i == 1 for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['#00FF00']),
                           marker="s", label=f"{vecteur}")
            elif all(i == 3 for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['white']),
                           marker="s", label=f"{vecteur}")
            else:
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['red', '#00FF00']),
                           marker="s", label=f"{vecteur}")
        else:
            if all(i == 0 for i in vecteur):
                ax.scatter(c_values, position * np.ones_like(c_values), c=vecteur, cmap=mcolors.ListedColormap(['red']),
                           marker="s", label=f"{vecteur}")
            elif all(i == 1 for i in vecteur):
                ax.scatter(c_values, position * np.ones_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['#00FF00']),
                           marker="s", label=f"{vecteur}")
            elif all(i == 3 for i in vecteur):
                ax.scatter(c_values, position * np.zeros_like(c_values), c=vecteur,
                           cmap=mcolors.ListedColormap(['white']),
                           marker="s", label=f"{vecteur}")
            else: ax.scatter(c_values, position * np.ones_like(c_values), c=vecteur, cmap=mcolors.ListedColormap(
                ['red', '#00FF00']),marker="s", label=f"{vecteur}")


def deviation_uni(symbole_initial, symbole_dev, vecteur_initiale, vecteur_deviation, vecteur_pays, structure, n,
                  c_values):
    global_var = {}
    for i in vecteur_pays:
        for x in range(len(symbole_dev)):
            global_var[f'w{i}_{symbole_initial}_{symbole_dev[x]}'] = []
    print(global_var)

    global_var_inverse = {}
    for i in vecteur_pays:
        for x in range(len(symbole_dev)):
            global_var_inverse[f'w{i}_{symbole_dev[x]}_{symbole_initial}'] = []
    print(global_var_inverse)

    for i in range(n):
        for x in range(len(vecteur_deviation)):
            for j in range(len(vecteur_pays)):
                key = f'w{vecteur_pays[j]}_{symbole_initial}_{symbole_dev[x]}'
                key_inverse = f'w{vecteur_pays[j]}_{symbole_dev[x]}_{symbole_initial}'
                if structure == 'bilatérale' and c_values[i] > 1 / 5:
                    global_var[key].append(2)
                elif (vecteur_initiale[j] - vecteur_deviation[x][j]).subs({c: c_values[i]}) > 0:
                    (global_var[key]).append(1)
                    (global_var_inverse[key_inverse]).append(0)
                else:
                    (global_var[key]).append(0)
                    (global_var_inverse[key_inverse]).append(1)

    global_dev = {}
    for x in range(len(symbole_dev)):
        global_dev[f'w_{symbole_initial}_{symbole_dev[x]}'] = []
    print(global_dev)

    global_dev_inverse = {}
    for x in range(len(symbole_dev)):
        global_dev_inverse[f'w_{symbole_dev[x]}_{symbole_initial}'] = []
    print(global_dev_inverse)

    for x in range(len(symbole_dev)):
        key = f'w_{symbole_initial}_{symbole_dev[x]}'
        global_dev[key].append(global_var[f'w{vecteur_pays[0]}_{symbole_initial}_{symbole_dev[x]}'])
        global_dev[key].append(global_var[f'w{vecteur_pays[1]}_{symbole_initial}_{symbole_dev[x]}'])
        global_dev[key].append(global_var[f'w{vecteur_pays[2]}_{symbole_initial}_{symbole_dev[x]}'])

    for x in range(len(symbole_dev)):
        key = f'w_{symbole_dev[x]}_{symbole_initial}'
        global_dev_inverse[key].append(global_var_inverse[f'w{vecteur_pays[0]}_{symbole_dev[x]}_{symbole_initial}'])
        global_dev_inverse[key].append(global_var_inverse[f'w{vecteur_pays[1]}_{symbole_dev[x]}_{symbole_initial}'])
        global_dev_inverse[key].append(global_var_inverse[f'w{vecteur_pays[2]}_{symbole_dev[x]}_{symbole_initial}'])


    print(global_var)
    print(global_var_inverse)

    return global_var, global_dev_inverse, global_dev


def deviation_coal(symbole_initial, symbole_dev, vecteur_initiale, vecteur_deviation, vecteur_credible,
                   vecteur_pays, structure, n, c_values):
    # créer les variables de déviation
    global_dev = {}
    for i in range(len(vecteur_deviation)):
        for j in vecteur_pays[i]:
            global_dev[f'w{j}_{symbole_initial}_{symbole_dev[i]}_{i+1}'] = []
    print(global_dev)


    for i in range(n):
        for x in range(len(vecteur_deviation)):
            for j in range(len(vecteur_pays[x])):
                key = f'w{vecteur_pays[x][j]}_{symbole_initial}_{symbole_dev[x]}_{x+1}'
                if structure == 'bilatérale' and c_values[i] > 1 / 5:
                    global_dev[key].append(2)
                elif len(vecteur_credible[x][j]) == 1 and vecteur_credible[x][j][0][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 2 and vecteur_credible[x][j][0][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 2 and vecteur_credible[x][j][1][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][0][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][1][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][2][i] == 1:
                    (global_dev[key]).append(1)
                elif vecteur_pays[x][j] == 'i' and (vecteur_initiale[0] - vecteur_deviation[x][0]).subs(
                        {c: c_values[i]}) > 0:
                    (global_dev[key]).append(1)
                elif vecteur_pays[x][j] == 'j' and (vecteur_initiale[1] - vecteur_deviation[x][1]).subs(
                        {c: c_values[i]}) > 0:
                    (global_dev[key]).append(1)
                elif vecteur_pays[x][j] == 'k' and (vecteur_initiale[2] - vecteur_deviation[x][2]).subs(
                        {c: c_values[i]}) > 0:
                    (global_dev[key]).append(1)
                else: (global_dev[key]).append(0)

    print(global_dev)

    return global_dev


##############################################################################################################
##############################################################################################################
########################### Fonction adapté pour reproduire Saggi ############################################
def deviation_uni_e(symbole_initial, symbole_dev, vecteur_initiale, vecteur_deviation, vecteur_pays, n):
    global_var = {}
    for i in vecteur_pays:
        for x in range(len(symbole_dev)):
            global_var[f'w{i}_{symbole_initial}_{symbole_dev[x]}'] = []
    print(global_var)

    global_var_inverse = {}
    for i in vecteur_pays:
        for x in range(len(symbole_dev)):
            global_var_inverse[f'w{i}_{symbole_dev[x]}_{symbole_initial}'] = []
    print(global_var_inverse)

    for i in range(n):
        for x in range(len(vecteur_deviation)):
            for j in range(len(vecteur_pays)):
                key = f'w{vecteur_pays[j]}_{symbole_initial}_{symbole_dev[x]}'
                key_inverse = f'w{vecteur_pays[j]}_{symbole_dev[x]}_{symbole_initial}'
                if (vecteur_initiale[j] - vecteur_deviation[x][j]).subs({e: 1}) > 0:
                    (global_var[key]).append(1)
                    (global_var_inverse[key_inverse]).append(0)
                else:
                    (global_var[key]).append(0)
                    (global_var_inverse[key_inverse]).append(1)

    global_dev = {}
    for x in range(len(symbole_dev)):
        global_dev[f'w_{symbole_initial}_{symbole_dev[x]}'] = []
    print(global_dev)

    global_dev_inverse = {}
    for x in range(len(symbole_dev)):
        global_dev_inverse[f'w_{symbole_dev[x]}_{symbole_initial}'] = []
    print(global_dev_inverse)

    for x in range(len(symbole_dev)):
        key = f'w_{symbole_initial}_{symbole_dev[x]}'
        global_dev[key].append(global_var[f'w{vecteur_pays[0]}_{symbole_initial}_{symbole_dev[x]}'])
        global_dev[key].append(global_var[f'w{vecteur_pays[1]}_{symbole_initial}_{symbole_dev[x]}'])
        global_dev[key].append(global_var[f'w{vecteur_pays[2]}_{symbole_initial}_{symbole_dev[x]}'])

    for x in range(len(symbole_dev)):
        key = f'w_{symbole_dev[x]}_{symbole_initial}'
        global_dev_inverse[key].append(global_var_inverse[f'w{vecteur_pays[0]}_{symbole_dev[x]}_{symbole_initial}'])
        global_dev_inverse[key].append(global_var_inverse[f'w{vecteur_pays[1]}_{symbole_dev[x]}_{symbole_initial}'])
        global_dev_inverse[key].append(global_var_inverse[f'w{vecteur_pays[2]}_{symbole_dev[x]}_{symbole_initial}'])


    print(global_var)
    print(global_var_inverse)

    return global_var, global_dev_inverse, global_dev
def deviation_coal_e(symbole_initial, symbole_dev, vecteur_initiale, vecteur_deviation, vecteur_credible,
                   vecteur_pays, n):
    # créer les variables de déviation
    global_dev = {}
    for i in range(len(vecteur_deviation)):
        for j in vecteur_pays[i]:
            global_dev[f'w{j}_{symbole_initial}_{symbole_dev[i]}_{i+1}'] = []
    print(global_dev)

    for i in range(n):
        for x in range(len(vecteur_deviation)):
            for j in range(len(vecteur_pays[x])):
                key = f'w{vecteur_pays[x][j]}_{symbole_initial}_{symbole_dev[x]}_{x + 1}'
                if len(vecteur_credible[x][j]) == 1 and vecteur_credible[x][j][0][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 2 and vecteur_credible[x][j][0][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 2 and vecteur_credible[x][j][1][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][0][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][1][i] == 1:
                    (global_dev[key]).append(1)
                elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][2][i] == 1:
                    (global_dev[key]).append(1)
                elif vecteur_pays[x] == 3:
                    if (vecteur_initiale[j] - vecteur_deviation[x][j]).subs({e: 1}) >= 0:
                        (global_dev[key]).append(1)
                elif vecteur_pays[x] == 2:
                    if (vecteur_initiale[j + 1] - vecteur_deviation[x][j + 1]).subs({e: 1}) > 0:
                        (global_dev[key]).append(1)
                else:
                    (global_dev[key]).append(0)

    print(global_dev)

    return global_dev

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
########################### Fonction adapté pour système d'équation linéaire ############################################
def obtenir_bornes_restrictives_3(s_t_o_b_i, s_t_o_b_i_dict, s_t_o_m_i, s_t_o_m_i_dict, X_b, X_m, variable_1,
                                  variable_2):
    def normaliser_inegalites(inegalites, variable):
        bornes = []
        for ineq in inegalites:
            if isinstance(ineq, Relational):  # Vérifier si c'est une inégalité
                if ineq.lhs == variable:  # Exemple: c <= 2 ou c >= -1
                    bornes.append(ineq) # append, car un seul élément
                elif ineq.rhs == variable:  # Exemple: -2/3 <= c  devient  c >= -2/3
                    bornes.append(ineq.reversed)
        return bornes
    exp_borne_b = []
    exp_borne_b.extend(s_t_o_b_i)
    exp_borne_b.extend(X_b)
    exp_borne_m = []
    exp_borne_m.extend(s_t_o_m_i)
    exp_borne_m.extend(X_m)
    result = s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_fij"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_fik"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_jh"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_kh"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_fij"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_fjk"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_ih"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_kh"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_fik"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_fjk"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_ih"]
    exp_borne_b.append(result)
    result = s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_jh"]
    exp_borne_b.append(result)
    result = s_t_o_m_i_dict["ti_quo"] - s_t_o_m_i_dict["ti_mij"]
    exp_borne_m.append(result)
    result = s_t_o_m_i_dict["ti_quo"] - s_t_o_m_i_dict["ti_mik"]
    exp_borne_m.append(result)
    result = s_t_o_m_i_dict["tj_quo"] - s_t_o_m_i_dict["tj_mij"]
    exp_borne_m.append(result)
    result = s_t_o_m_i_dict["tj_quo"] - s_t_o_m_i_dict["tj_mjk"]
    exp_borne_m.append(result)
    result = s_t_o_m_i_dict["tk_quo"] - s_t_o_m_i_dict["tk_mik"]
    exp_borne_m.append(result)
    result = s_t_o_m_i_dict["tk_quo"] - s_t_o_m_i_dict["tk_mjk"]
    exp_borne_m.append(result)
    bornes_bilaterales_ci = []
    bornes_multilaterales_ci = []
    bornes_bilaterales_cj = []
    bornes_multilaterales_cj = []
    bornes_b_ci_max = []
    bornes_b_ci_min = []
    bornes_b_cj_max = []
    bornes_b_cj_min = []
    bornes_m_ci_max = []
    bornes_m_ci_min = []
    bornes_m_cj_max = []
    bornes_m_cj_min = []
    bornes_restrictives_bilaterale_ci = []
    bornes_restrictives_multilaterale_ci = []
    bornes_restrictives_bilaterale_cj = []
    bornes_restrictives_multilaterale_cj = []
    for tarif in range(len(exp_borne_b)):
        if exp_borne_b[tarif].free_symbols == {ci}:
            result = solve(exp_borne_b[tarif]>=0, variable_1)
            bornes_bilaterales_ci.extend(result.args) # .extend() , car solve() renvoit 2 conditions en même temps
    for tarif in range(len(exp_borne_m)):
        if exp_borne_b[tarif].free_symbols == {ci}:
            result = solve(exp_borne_m[tarif]>=0, variable_1)
            bornes_multilaterales_ci.extend(result.args)
    bornes_bilaterales_ci = normaliser_inegalites(bornes_bilaterales_ci, ci)
    bornes_multilaterales_ci = normaliser_inegalites(bornes_multilaterales_ci, ci)

    for borne in bornes_bilaterales_ci:
        if isinstance(borne, Relational) and borne.rel_op == '<=':
            valeur = borne.rhs
            bornes_b_ci_max.append(valeur)

    for borne in bornes_bilaterales_ci:
        if isinstance(borne, Relational) and borne.rel_op == '>=':
            valeur = borne.rhs
            bornes_b_ci_min.append(valeur)
    bornes_restrictives_bilaterale_ci.append(Max(*bornes_b_ci_min) if Max(*bornes_b_ci_min) > 0 else 0)
    bornes_restrictives_bilaterale_ci.append(Min(*bornes_b_ci_max))

    for borne in bornes_multilaterales_ci:
        if isinstance(borne, Relational) and borne.rel_op == '<=':
            valeur = borne.rhs
            bornes_m_ci_max.append(valeur)

    for borne in bornes_multilaterales_ci:
        if isinstance(borne, Relational) and borne.rel_op == '>=':
            valeur = borne.rhs
            bornes_m_ci_min.append(valeur)

    bornes_restrictives_multilaterale_ci.append(Max(*bornes_m_ci_min) if Max(*bornes_m_ci_min) > 0 else 0)
    bornes_restrictives_multilaterale_ci.append(Min(*bornes_m_ci_max))
    print(f'borne restrictive bilatérale ci = {bornes_restrictives_bilaterale_ci}\n')
    print(f'borne restrictive multilatérale ci = {bornes_restrictives_multilaterale_ci}\n')
    
    for tarif in range(len(exp_borne_b)):
        if exp_borne_b[tarif].free_symbols == {cj}:
            result = solve(exp_borne_b[tarif]>=0, variable_2)
            bornes_bilaterales_cj.extend(result.args) # .extend() , car solve() renvoit 2 conditions en même temps
    for tarif in range(len(exp_borne_m)):
        if exp_borne_b[tarif].free_symbols == {cj}:
            result = solve(exp_borne_m[tarif]>=0, variable_2)
            bornes_multilaterales_cj.extend(result.args)
    bornes_bilaterales_cj = normaliser_inegalites(bornes_bilaterales_cj, cj)
    bornes_multilaterales_cj = normaliser_inegalites(bornes_multilaterales_cj, cj)

    for borne in bornes_bilaterales_cj:
        if isinstance(borne, Relational) and borne.rel_op == '<=':
            valeur = borne.rhs
            bornes_b_cj_max.append(valeur)

    for borne in bornes_bilaterales_cj:
        if isinstance(borne, Relational) and borne.rel_op == '>=':
            valeur = borne.rhs
            bornes_b_cj_min.append(valeur)
    bornes_restrictives_bilaterale_cj.append(Max(*bornes_b_cj_min) if Max(*bornes_b_cj_min) > 0 else 0)
    bornes_restrictives_bilaterale_cj.append(Min(*bornes_b_cj_max))

    for borne in bornes_multilaterales_cj:
        if isinstance(borne, Relational) and borne.rel_op == '<=':
            valeur = borne.rhs
            bornes_m_cj_max.append(valeur)

    for borne in bornes_multilaterales_cj:
        if isinstance(borne, Relational) and borne.rel_op == '>=':
            valeur = borne.rhs
            bornes_m_cj_min.append(valeur)

    bornes_restrictives_multilaterale_cj.append(Max(*bornes_m_cj_min) if Max(*bornes_m_cj_min) > 0 else 0)
    bornes_restrictives_multilaterale_cj.append(Min(*bornes_m_cj_max))
    print(f'borne restrictive bilatérale cj = {bornes_restrictives_bilaterale_cj}\n')
    print(f'borne restrictive multilatérale cj = {bornes_restrictives_multilaterale_cj}\n')

def deviation_uni_matrice(symbole_initial, symbole_dev, vecteur_initiale, vecteur_deviation, vecteur_pays, structure,
                          n_1, n_2, ci_values, cj_values):
    global_var = {}
    for i in vecteur_pays:
        for x in range(len(symbole_dev)):
            global_var[f'w{i}_{symbole_initial}_{symbole_dev[x]}'] = np.zeros((n_1, n_2))
    print(global_var)

    global_var_inverse = {}
    for i in vecteur_pays:
        for x in range(len(symbole_dev)):
            global_var_inverse[f'w{i}_{symbole_dev[x]}_{symbole_initial}'] = np.zeros((n_1, n_2))
    print(global_var_inverse)

    for i in range(n_1):
        for j in range(n_2):
            for x in range(len(vecteur_deviation)):
                for y in range(len(vecteur_pays)):
                    key = f'w{vecteur_pays[y]}_{symbole_initial}_{symbole_dev[x]}'
                    key_inverse = f'w{vecteur_pays[y]}_{symbole_dev[x]}_{symbole_initial}'
                    if structure == 'bilatérale' and cj_values[j] > 1/5:
                        global_var[key][i, j] = 2
                    elif (vecteur_initiale[y] - vecteur_deviation[x][y]).subs({ci: ci_values[i], cj: cj_values[j]}) > 0:
                        global_var[key][i, j] = 1
                        global_var_inverse[key_inverse][i, j] = 0
                    else:
                        global_var[key][i, j] = 0
                        global_var_inverse[key_inverse][i, j] = 1


    print(global_var)
    print(global_var_inverse)

    return global_var, global_var_inverse


def deviation_coal_matrice(symbole_initial, symbole_dev, vecteur_initiale, vecteur_deviation, vecteur_credible,
                   vecteur_pays, structure, n_1, n_2, ci_values, cj_values):
    # créer les variables de déviation
    global_dev = {}
    for i in range(len(vecteur_deviation)):
        for j in vecteur_pays[i]:
            global_dev[f'w{j}_{symbole_initial}_{symbole_dev[i]}_{i+1}'] = np.zeros((n_1, n_2))
    print(global_dev)

    for i in range(n_1):
        for k in range(n_2):
            for x in range(len(vecteur_deviation)):
                for j in range(len(vecteur_pays[x])):
                    key = f'w{vecteur_pays[x][j]}_{symbole_initial}_{symbole_dev[x]}_{x + 1}'
                    if structure == 'bilatérale' and cj_values[k] > 1/5:
                        global_dev[key][i, k] = 2
                    elif len(vecteur_credible[x][j]) == 1 and vecteur_credible[x][j][0][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 2 and vecteur_credible[x][j][0][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 2 and vecteur_credible[x][j][1][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][0][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][1][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][2][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif vecteur_pays[x][j] == 'i' and (vecteur_initiale[0] - vecteur_deviation[x][0]).subs(
                            {ci: ci_values[i], cj: cj_values[k]}) > 0:
                        global_dev[key][i, k] = 1
                    elif vecteur_pays[x][j] == 'j' and (vecteur_initiale[1] - vecteur_deviation[x][1]).subs(
                            {ci: ci_values[i], cj: cj_values[k]}) > 0:
                        global_dev[key][i, k] = 1
                    elif vecteur_pays[x][j] == 'k' and (vecteur_initiale[2] - vecteur_deviation[x][2]).subs(
                            {ci: ci_values[i], cj: cj_values[k]}) > 0:
                        global_dev[key][i, k] = 1
                    else:
                        global_dev[key][i, k] = 0


    print(global_dev)

    return global_dev

def graph_matrice(structure, n_b_ci, n_b_cj, matrice_dev, titre):
    tick_positions_1 = np.linspace(0, n_b_ci, num=5)
    tick_positions_2 = np.linspace(0, n_b_cj, num=5)
    tick_etiquette_1 = np.round(np.linspace(0, 1, num=5), 2)
    tick_etiquette_2 = np.round(np.linspace(0, 2/5, num=5), 2)
    plt.xticks(tick_positions_2, tick_etiquette_2)
    plt.yticks(tick_positions_1, tick_etiquette_1)
    # Vérification des conditions sur la matrice
    if structure == "bilatérale":
        unique_values = set(np.unique(matrice_dev))  # Trouve les valeurs uniques présentes dans la matrice
        if unique_values == {0}:
            cmap = mcolors.ListedColormap(['red'])
        elif unique_values == {1}:
            cmap = mcolors.ListedColormap(['#00FF00'])
        elif unique_values == {2}:
            cmap = mcolors.ListedColormap(['black'])
        elif unique_values == {0, 1}:
            cmap = mcolors.ListedColormap(['red', '#00FF00'])
        elif unique_values == {0, 2}:
            cmap = mcolors.ListedColormap(['red', 'black'])
        elif unique_values == {1, 2}:
            cmap = mcolors.ListedColormap(['#00FF00', 'black'])
        elif unique_values == {0, 1, 2}:
            cmap = mcolors.ListedColormap(['red', '#00FF00', 'black'])
    else:
        if np.all(matrice_dev == 0):
            cmap = mcolors.ListedColormap(['red'])
        elif np.all(matrice_dev == 1):
            cmap = mcolors.ListedColormap(['#00FF00'])
        else:
            cmap = mcolors.ListedColormap(['red', '#00FF00'])



    # Affichage de la matrice
    plt.imshow(matrice_dev, cmap=cmap, origin="lower")
    plt.xlabel("Coûts de tansport $\lambda_{l}$")
    plt.ylabel("Coûts de tansport $\lambda_{i}$")
    plt.title(titre)
    plt.show()

def graph_matrice_e(structure, n_b_ci, n_b_cj, titre, situation_initiale, matrice_dict):
    # Vérification des conditions sur la matrice
    for nom, matrice in matrice_dict.items():
        if structure == "bilatérale":
            unique_values = set(np.unique(matrice))  # Trouve les valeurs uniques présentes dans la matrice
            if unique_values == {0}:
                cmap = mcolors.ListedColormap(['red'])
            elif unique_values == {1}:
                cmap = mcolors.ListedColormap(['#00FF00'])
            elif unique_values == {2}:
                cmap = mcolors.ListedColormap(['black'])
            elif unique_values == {0, 1}:
                cmap = mcolors.ListedColormap(['red', '#00FF00'])
            elif unique_values == {0, 2}:
                cmap = mcolors.ListedColormap(['red', 'black'])
            elif unique_values == {1, 2}:
                cmap = mcolors.ListedColormap(['#00FF00', 'black'])
            elif unique_values == {0, 1, 2}:
                cmap = mcolors.ListedColormap(['red', '#00FF00', 'black'])
        else:
            if np.all(matrice == 0):
                cmap = mcolors.ListedColormap(['red'])
            elif np.all(matrice == 1):
                cmap = mcolors.ListedColormap(['#00FF00'])
            else:
                cmap = mcolors.ListedColormap(['red', '#00FF00'])
        tick_positions_1 = np.linspace(0, n_b_ci, num=5)
        tick_positions_2 = np.linspace(0, n_b_cj, num=5)
        tick_etiquette_1 = np.round(np.linspace(0, 1, num=5), 2)
        tick_etiquette_2 = np.round(np.linspace(0, 2/5, num=5), 2)
        plt.xticks(tick_positions_2, tick_etiquette_2)
        plt.yticks(tick_positions_1, tick_etiquette_1)
        plt.imshow(matrice, cmap=cmap, origin="lower")
        plt.xlabel(r"Coûts de tansport $\lambda_{l}$")
        plt.ylabel(r"Coûts de tansport $\lambda_{i}$")
        plt.title(f'déviation {titre} {situation_initiale} - {nom}')
        plt.show()

def graph_matrice_presentation(n_b_ci, n_b_cj, matrice_dev, titre):
    tick_positions_1 = np.linspace(0, n_b_ci, num=5)
    tick_positions_2 = np.linspace(0, n_b_cj, num=5)
    tick_etiquette_1 = np.round(np.linspace(0, 1, num=5), 2)
    tick_etiquette_2 = np.round(np.linspace(0, 2/5, num=5), 2)
    plt.xticks(tick_positions_1, tick_etiquette_1)
    plt.yticks(tick_positions_2, tick_etiquette_2)

    cmap = mcolors.ListedColormap(['red', '#00FF00', 'cyan', 'magenta'])

    # Affichage de la matrice
    plt.imshow(matrice_dev, cmap=cmap, origin="lower")
    plt.xlabel("Coûts de tansport $\lambda_{i}$")
    plt.ylabel("Coûts de tansport $\lambda_{l}$")
    plt.title(titre)

    legend_elements = [
        Patch(facecolor='#00FF00', edgecolor='black', label='ENEC_$F^m$'),
        Patch(facecolor='cyan', edgecolor='black', label='ENEC_F'),
        Patch(facecolor='magenta', edgecolor='black', label='F et $F^m$ ENEC'),
        Patch(facecolor='red', edgecolor='black', label='Aucun ENEC')
    ]
    plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1))

    plt.tight_layout()

    plt.show()


##############################################################################################################
##############################################################################################################
##############################################################################################################
################################### situation 4 #########################################

def mettre_en_fonction_de_c(s_t_o_b_i_dict, s_t_o_m_i_dict, s_t_o_b_i_vect, s_t_o_m_i_vect, X_b, X_m, p_values):
    s_t_o_b_i_all = []
    s_t_o_m_i_all = []
    X_b_all = []
    X_m_all = []
    tarif_b_all = []
    tarif_m_all = []
    tarif_b_all = []
    tarif_m_all = []
    conditions_T_m_all = []
    conditions_T_b_all = []

    result = (s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_fij"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_fik"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_jh"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["ti_quo"] - s_t_o_b_i_dict["ti_kh"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_fij"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_fjk"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_ih"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["tj_quo"] - s_t_o_b_i_dict["tj_kh"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_fik"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_fjk"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_ih"])
    conditions_T_b_all.append(result)
    result = (s_t_o_b_i_dict["tk_quo"] - s_t_o_b_i_dict["tk_jh"])
    conditions_T_b_all.append(result)

    result = (s_t_o_m_i_dict["ti_quo"] - s_t_o_m_i_dict["ti_mij"])
    conditions_T_m_all.append(result)
    result = (s_t_o_m_i_dict["ti_quo"] - s_t_o_m_i_dict["ti_mik"])
    conditions_T_m_all.append(result)
    result = (s_t_o_m_i_dict["tj_quo"] - s_t_o_m_i_dict["tj_mij"])
    conditions_T_m_all.append(result)
    result = (s_t_o_m_i_dict["tj_quo"] - s_t_o_m_i_dict["tj_mjk"])
    conditions_T_m_all.append(result)
    result = (s_t_o_m_i_dict["tk_quo"] - s_t_o_m_i_dict["tk_mik"])
    conditions_T_m_all.append(result)
    result = (s_t_o_m_i_dict["tk_quo"] - s_t_o_m_i_dict["tk_mjk"])
    conditions_T_m_all.append(result)

    for k in range(len(p_values)):
        tarif_b_all.append([eq.subs({p: p_values[k]}) for eq in conditions_T_b_all])
        s_t_o_b_i_all.append([eq.subs({p: p_values[k]}) for eq in s_t_o_b_i_vect])
        X_b_all.append([eq.subs({p: p_values[k]}) for eq in X_b])
        tarif_m_all.append([eq.subs({p: p_values[k]}) for eq in conditions_T_m_all])
        s_t_o_m_i_all.append([eq.subs({p: p_values[k]}) for eq in s_t_o_m_i_vect])
        X_m_all.append([eq.subs({p: p_values[k]}) for eq in X_m])

    return (s_t_o_b_i_all, X_b_all, tarif_b_all, s_t_o_m_i_all, X_m_all, tarif_m_all)


def trier_les_equations_et_solver(s_t_o_b_i_all, X_b_all, tarif_b_all, s_t_o_m_i_all, X_m_all, tarif_m_all):
    liste_solutions_bilaterales_c = []
    liste_solutions_multilaterales_c = []
    liste_eq_sans_solution_b = []
    liste_eq_sans_solution_m = []
    liste_problematique_b = []
    liste_problematique_m = []

    for eq_list in s_t_o_b_i_all:
        solution = []
        for eq in eq_list:
            try:
                solution = simplify(sp.solve(eq >= 0, c))
                if not solution:
                    liste_eq_sans_solution_b.append(eq)
                else:
                    liste_solutions_bilaterales_c.append(solution)
            except Exception as e:
                liste_problematique_b.append(eq)

    for eq_list in X_b_all:
        solution = []
        for eq in eq_list:
            try:
                solution = simplify(sp.solve(eq >= 0, c))
                if not solution:
                    liste_eq_sans_solution_b.append(eq)
                else:
                    liste_solutions_bilaterales_c.append(solution)
            except Exception as e:
                liste_problematique_b.append(eq)

    for eq_list in tarif_b_all:
        solution = []
        for eq in eq_list:
            try:
                solution = simplify(sp.solve(eq >= 0, c))
                if not solution:
                    liste_eq_sans_solution_b.append(eq)
                else:
                    liste_solutions_bilaterales_c.append(solution)
            except Exception as e:
                liste_problematique_b.append(eq)

    for eq_list in s_t_o_m_i_all:
        solution = []
        for eq in eq_list:
            try:
                solution = simplify(sp.solve(eq >= 0, c))
                if not solution:
                    liste_eq_sans_solution_m.append(eq)
                else:
                    liste_solutions_multilaterales_c.append(solution)
            except Exception as e:
                liste_problematique_m.append(eq)

    for eq_list in X_m_all:
        solution = []
        for eq in eq_list:
            try:
                solution = simplify(sp.solve(eq >= 0, c))
                if not solution:
                    liste_eq_sans_solution_m.append(eq)
                else:
                    liste_solutions_multilaterales_c.append(solution)
            except Exception as e:
                liste_problematique_m.append(eq)

    for eq_list in tarif_m_all:
        solution = []
        for eq in eq_list:
            try:
                solution = simplify(sp.solve(eq >= 0, c))
                if not solution:
                    liste_eq_sans_solution_m.append(eq)
                else:
                    liste_solutions_multilaterales_c.append(solution)
            except Exception as e:
                liste_problematique_m.append(eq)

    for eq in liste_problematique_b[:]:
        print(type(eq), eq)

        solutions = []
        try:
            ineq = eq >= 0
            solution = reduce_inequalities(ineq)

            if solution:
                solutions.append(solution)
                liste_problematique_b.remove(eq)

        except Exception as e:
            pass

    for eq in liste_problematique_m[:]:  # Copier la liste pour éviter les erreurs de modification en place
        print(type(eq), eq)  # Debugging

        solutions = []
        try:
            ineq = eq >= 0  # Construire explicitement l'inégalité
            solution = reduce_inequalities(ineq)

            if solution:  # Si une solution a été trouvée
                solutions.append(solution)
                liste_problematique_m.remove(eq)  # Supprime l'équation de la liste

        except Exception as e:
            pass

    print("Longueurs des vecteurs problematiques pour vérification ")
    print(len(liste_problematique_b))
    print(len(liste_problematique_m))

    return (liste_solutions_bilaterales_c, liste_problematique_b, liste_eq_sans_solution_b,
            liste_solutions_multilaterales_c, liste_problematique_m, liste_eq_sans_solution_m)


def normaliser_inegalites(inegalites):
    bornes = []
    for ineq in inegalites:
        if isinstance(ineq, Relational):  # Vérifier si c'est une inégalité
            if ineq.lhs == c:
                bornes.append(ineq)  # Ajouter telle quelle
            elif ineq.rhs == c:
                bornes.append(ineq.reversed)  # Inverser l'inégalité
        elif isinstance(ineq, sp.And):  # Si l'inégalité est une conjonction (et)
            # Séparer les inégalités dans la conjonction
            bornes.extend(
                normaliser_inegalites(list(ineq.args)))  # On utilise .args pour itérer sur les éléments de l'And
    return bornes


def extraire_bornes_bilaterales(inegalites_normalisees):
    bornes_b_min = []
    bornes_b_max = []
    bornes_restrictives_bilaterale = []

    for ineq in inegalites_normalisees:
        if isinstance(ineq, Relational) and ineq.rel_op == '<=':
            valeur = ineq.rhs
            bornes_b_max.append(valeur)

    for ineq in inegalites_normalisees:
        if isinstance(ineq, Relational) and ineq.rel_op == '>=':
            valeur = ineq.rhs
            bornes_b_min.append(valeur)

    bornes_restrictives_bilaterale.append(Max(*bornes_b_min) if Max(*bornes_b_min) > 0 else 0)
    bornes_restrictives_bilaterale.append(Min(*bornes_b_max))

    print(f'borne restrictive bilatérale = {bornes_restrictives_bilaterale}\n')
    return bornes_restrictives_bilaterale


def extraire_bornes_multilaterales(inegalites_normalisees):
    bornes_m_min = []
    bornes_m_max = []
    bornes_restrictives_multilaterale = []

    for ineq in inegalites_normalisees:
        if isinstance(ineq, Relational) and ineq.rel_op == '<=':
            valeur = ineq.rhs
            bornes_m_max.append(valeur)

    for ineq in inegalites_normalisees:
        if isinstance(ineq, Relational) and ineq.rel_op == '>=':
            valeur = ineq.rhs
            bornes_m_min.append(valeur)

    bornes_restrictives_multilaterale.append(Max(*bornes_m_min) if Max(*bornes_m_min) > 0 else 0)
    bornes_restrictives_multilaterale.append(Min(*bornes_m_max))

    print(f'borne restrictive multilatérale = {bornes_restrictives_multilaterale}\n')
    return bornes_restrictives_multilaterale


def deviation_uni_matrice_pythagore(symbole_initial, symbole_dev, vecteur_initiale, vecteur_deviation, vecteur_pays,
                                    structure, n_1, n_2, ci_values, p_values):
    global_var = {}
    for i in vecteur_pays:
        for x in range(len(symbole_dev)):
            global_var[f'w{i}_{symbole_initial}_{symbole_dev[x]}'] = np.zeros((n_1, n_2))
    print(global_var)

    global_var_inverse = {}
    for i in vecteur_pays:
        for x in range(len(symbole_dev)):
            global_var_inverse[f'w{i}_{symbole_dev[x]}_{symbole_initial}'] = np.zeros((n_1, n_2))
    print(global_var_inverse)

    for i in range(n_1):
        for j in range(n_2):
            for x in range(len(vecteur_deviation)):
                for y in range(len(vecteur_pays)):
                    key = f'w{vecteur_pays[y]}_{symbole_initial}_{symbole_dev[x]}'
                    key_inverse = f'w{vecteur_pays[y]}_{symbole_dev[x]}_{symbole_initial}'
                    if structure == 'bilatérale' and ci_values[i] > 1 / 5:
                        global_var[key][i, j] = 2
                    elif (vecteur_initiale[y] - vecteur_deviation[x][y]).subs({c: ci_values[i], p: p_values[j]}) >= 0:
                        global_var[key][i, j] = 1
                        global_var_inverse[key_inverse][i, j] = 0
                    else:
                        global_var[key][i, j] = 0
                        global_var_inverse[key_inverse][i, j] = 1

    print(global_var)
    print(global_var_inverse)

    return global_var, global_var_inverse


def deviation_coal_matrice_pythagore(symbole_initial, symbole_dev, vecteur_initiale, vecteur_deviation,
                                     vecteur_credible, vecteur_pays, structure, n_1, n_2, ci_values, p_values):
    # créer les variables de déviation
    global_dev = {}
    for i in range(len(vecteur_deviation)):
        for j in vecteur_pays[i]:
            global_dev[f'w{j}_{symbole_initial}_{symbole_dev[i]}_{i + 1}'] = np.zeros((n_1, n_2))
    print(global_dev)

    for i in range(n_1):
        for k in range(n_2):
            for x in range(len(vecteur_deviation)):
                for j in range(len(vecteur_pays[x])):
                    key = f'w{vecteur_pays[x][j]}_{symbole_initial}_{symbole_dev[x]}_{x + 1}'
                    if structure == 'bilatérale' and ci_values[i] > 1 / 5:
                        global_dev[key][i, k] = 2
                    elif len(vecteur_credible[x][j]) == 1 and vecteur_credible[x][j][0][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 2 and vecteur_credible[x][j][0][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 2 and vecteur_credible[x][j][1][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][0][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][1][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif len(vecteur_credible[x][j]) == 3 and vecteur_credible[x][j][2][i, k] == 1:
                        global_dev[key][i, k] = 1
                    elif vecteur_pays[x][j] == 'i' and (vecteur_initiale[0] - vecteur_deviation[x][0]).subs(
                            {c: ci_values[i], p: p_values[k]}) > 0:
                        global_dev[key][i, k] = 1
                    elif vecteur_pays[x][j] == 'j' and (vecteur_initiale[1] - vecteur_deviation[x][1]).subs(
                            {c: ci_values[i], p: p_values[k]}) > 0:
                        global_dev[key][i, k] = 1
                    elif vecteur_pays[x][j] == 'k' and (vecteur_initiale[2] - vecteur_deviation[x][2]).subs(
                            {c: ci_values[i], p: p_values[k]}) > 0:
                        global_dev[key][i, k] = 1
                    else:
                        global_dev[key][i, k] = 0

    print(global_dev)

    return global_dev


def graph_matrice_e_pythagore(n_b_cij, n_p, titre, situation_initiale, matrice_dict):
    for nom, matrice in matrice_dict.items():
        # Définition de la colormap en fonction du contenu de la matrice
        if np.all(matrice == 0):
            cmap = mcolors.ListedColormap(['red'])
        elif np.all(matrice == 1):
            cmap = mcolors.ListedColormap(['#00FF00'])
        elif np.all(matrice == 2):
            cmap = mcolors.ListedColormap(['black'])
        elif np.all(np.isin(matrice, [0, 1])):
            cmap = mcolors.ListedColormap(['red', '#00FF00'])
        elif np.all(np.isin(matrice, [1, 2])):
            cmap = mcolors.ListedColormap(['#00FF00', 'black'])
        elif np.all(np.isin(matrice, [0, 2])):
            cmap = mcolors.ListedColormap(['red', 'black'])
        else:
            cmap = mcolors.ListedColormap(['red', '#00FF00', 'black'])

        # Création des ticks
        tick_positions_1 = np.linspace(0, n_p, num=5)
        tick_positions_2 = np.linspace(0, n_b_cij, num=5)
        tick_etiquette_1 = np.round(np.linspace(0, 1, num=5), 2)
        tick_etiquette_2 = np.round(np.linspace(0, 2 / 5, num=5), 2)

        # Création de la figure
        plt.figure()
        plt.xticks(tick_positions_1, tick_etiquette_1)
        plt.yticks(tick_positions_2, tick_etiquette_2)
        plt.imshow(matrice, cmap=cmap, origin="lower")
        plt.xlabel("Proportion (p)")
        plt.ylabel(r'$\lambda_{ij}$')
        plt.title(f'déviation {titre} {situation_initiale} - {nom}')

        # Affichage de la figure
        plt.show()


def graph_matrice_pythagore(n_b_cij, n_p, matrice_dev, titre):
    tick_positions_1 = np.linspace(0, n_p, num=5)
    tick_positions_2 = np.linspace(0, n_b_cij, num=5)
    tick_etiquette_1 = np.round(np.linspace(0, 1, num=5), 2)
    tick_etiquette_2 = np.round(np.linspace(0, 2 / 5, num=5), 2)
    plt.xticks(tick_positions_1, tick_etiquette_1)
    plt.yticks(tick_positions_2, tick_etiquette_2)
    # Vérification des conditions sur la matrice
    if np.all(matrice_dev == 0):
        cmap = mcolors.ListedColormap(['red'])
    elif np.all(matrice_dev == 1):
        cmap = mcolors.ListedColormap(['#00FF00'])
    elif np.all(matrice_dev == 2):
        cmap = mcolors.ListedColormap(['black'])
    elif np.all(np.isin(matrice_dev, [0, 1])):
        cmap = mcolors.ListedColormap(['red', '#00FF00'])
    elif np.all(np.isin(matrice_dev, [1, 2])):
        cmap = mcolors.ListedColormap(['#00FF00', 'black'])
    elif np.all(np.isin(matrice_dev, [0, 2])):
        cmap = mcolors.ListedColormap(['red', 'black'])
    else:
        cmap = mcolors.ListedColormap(['red', '#00FF00', 'black'])

    # Affichage de la matrice
    plt.imshow(matrice_dev, cmap=cmap, origin="lower")
    plt.xlabel("Proportion (p)")
    plt.ylabel(r'$\lambda_{ij}$')
    plt.title(titre)

    plt.show()