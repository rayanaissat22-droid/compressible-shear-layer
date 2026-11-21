import os
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import yaml

DOSSIER = "Outputs"


class CoucheCisaillement:

    def __init__(self, lam, Pr, T1, T2, L_haut, L_zero, h):
        self.lam = lam
        self.Pr = Pr
        self.T1 = T1
        self.T2 = T2
        self.h = h
        self.L_haut = L_haut
        if lam > 0:
            self.L_bas = np.ceil(L_haut / np.sqrt(lam))   # arrondi à l'entier supérieur
        else:
            self.L_bas = L_zero

    def systeme_F(self, x):
        f, u, g = x
        return np.array([u, g, -f * g])

    def systeme_FG(self, x):
        f, u, g, H, dH = x
        return np.array([u, g, -f * g, dH, -self.Pr * f * dH])

    def pas_rk4(self, systeme, x, h):
        k1 = systeme(x)
        k2 = systeme(x + h/2 * k1)
        k3 = systeme(x + h/2 * k2)
        k4 = systeme(x + h * k3)
        xn = x + (h/6) * (k1 + 2 * k2 + 2 * k3 + k4)
        return xn

    def integrer(self, systeme, x0, eta_fin):
        h = self.h if eta_fin > 0 else -self.h
        n = int(round(eta_fin / h))
        eta = np.linspace(0, eta_fin, n + 1)
        X = np.zeros((n + 1, len(x0)))
        X[0] = x0
        for i in range(n):
            X[i + 1] = self.pas_rk4(systeme, X[i], h)
        return eta, X

    def residus(self, s):
        eta, X = self.integrer(self.systeme_F, s, self.L_haut)
        eta2, X2 = self.integrer(self.systeme_F, s, -self.L_bas)
        f_haut, u_haut, g_haut = X[-1]
        f_bas, u_bas, g_bas = X2[-1]
        D1 = u_haut - 1
        D2 = u_bas - self.lam
        D3 = (self.L_haut*u_haut - f_haut) + self.lam * (-self.L_bas*u_bas - f_bas)
        return np.array([D1, D2, D3])

    def jacobienne(self, s, eps):
        D0 = self.residus(s)
        J = np.zeros((3, 3))
        for j in range(3):
            s_bouge = s.copy()
            s_bouge[j] = s_bouge[j] + eps
            D_bouge = self.residus(s_bouge)
            J[:, j] = (D_bouge - D0) / eps
        return J

    def newton(self, s0, eps, tol, itmax):
        s = s0.copy()
        for it in range(itmax):
            D = self.residus(s)
            norme = np.linalg.norm(D)
            if norme < tol:
                return s
            J = self.jacobienne(s, eps)
            ds = np.linalg.solve(J, -D)
            s = s + ds
        print("Attention : Newton n'a pas convergé")
        return s

    def profil_complet(self, s):
        eta_h, X_h = self.integrer(self.systeme_F, s, self.L_haut)
        eta_b, X_b = self.integrer(self.systeme_F, s, -self.L_bas)
        eta = np.concatenate((eta_b[::-1], eta_h[1:]))
        X = np.concatenate((X_b[::-1], X_h[1:]))
        return eta, X

    def resoudre_F(self, s0, eps, tol, itmax):
        self.s = self.newton(s0, eps, tol, itmax)
        self.eta, self.X = self.profil_complet(self.s)

    def resoudre_G(self):
        f0, u0, g0 = self.s
        x_tir1 = np.array([f0, u0, g0, 1, 0])
        x_tir2 = np.array([f0, u0, g0, 0, 1])

        eta_h, X1_h = self.integrer(self.systeme_FG, x_tir1, self.L_haut)
        eta_b, X1_b = self.integrer(self.systeme_FG, x_tir1, -self.L_bas)
        eta_h, X2_h = self.integrer(self.systeme_FG, x_tir2, self.L_haut)
        eta_b, X2_b = self.integrer(self.systeme_FG, x_tir2, -self.L_bas)

        M = np.array([[X1_b[-1, 3], X2_b[-1, 3]],
                      [X1_h[-1, 3], X2_h[-1, 3]]])
        self.a, self.b = np.linalg.solve(M, [0, 1])

        X1 = np.concatenate([X1_b[::-1], X1_h[1:]])
        X2 = np.concatenate([X2_b[::-1], X2_h[1:]])
        self.H = self.a * X1[:, 3] + self.b * X2[:, 3]       # G'
        self.dH = self.a * X1[:, 4] + self.b * X2[:, 4]      # G''
        self.G = integrale(self.H, self.eta)                 # G, avec G(0) = 0

    def ligne_partage(self):
        F = self.X[:, 0]
        for i in range(len(F) - 1):
            if F[i] < 0 and F[i + 1] >= 0:
                eta_i = np.interp(0, [F[i], F[i + 1]], [self.eta[i], self.eta[i + 1]])
                u_i = np.interp(eta_i, self.eta, self.X[:, 1])
                return eta_i, u_i

    def epaisseur_vorticite(self):
        return (1 - self.lam) / np.max(self.X[:, 2])

    def vitesses_infini(self):
        V_haut = self.eta[-1] * self.X[-1, 1] - self.X[-1, 0]
        V_bas = self.eta[0] * self.X[0, 1] - self.X[0, 0]
        return V_haut / np.sqrt(2), V_bas / np.sqrt(2)

    def compressible(self):
        # température et masse volumique réduite
        self.T = self.T2 + (self.T1 - self.T2) * self.H
        self.rho_s = self.T1 / self.T

        # coordonnée physique y/delta (Howarth inversé), avec Y = 0 en eta = 0
        self.Y = integrale(self.T / self.T1, self.eta)

        # grandeurs compressibles
        self.rho0 = np.interp(0, self.eta, self.rho_s)
        self.delta_w = (1 - self.lam) / np.max(self.rho_s * self.X[:, 2])
        self.Cf_c = self.s[2] * self.delta_w
        self.Nu_c = -self.rho0 * self.b * self.delta_w

    def calculer(self, s0, eps, tol, itmax):
        self.resoudre_F(s0, eps, tol, itmax)
        self.resoudre_G()

        # 2) grandeurs en variable eta
        eta_i, u_i = self.ligne_partage()
        eta_omega = self.epaisseur_vorticite()
        V_pinf, V_minf = self.vitesses_infini()
        Cf = self.s[2] * eta_omega
        Nu = -self.b * eta_omega

        # 3) retour au compressible
        self.compressible()

        # 4) les deux lignes de résultats
        self.ligne_eta = [self.lam, self.s[0], self.s[1], self.s[2], eta_i, u_i, eta_omega,
                          V_pinf, V_minf, Cf, self.a, self.b, Nu, -self.L_bas, self.L_haut]
        self.ligne_comp = [self.lam, self.rho0, self.delta_w, self.Cf_c, self.Nu_c]

        return self.s

    def tracer_F(self):
        plt.figure()
        plt.plot(self.eta, self.X[:, 0], label="F")
        plt.plot(self.eta, self.X[:, 1], label="F'")
        plt.plot(self.eta, self.X[:, 2], label="F''")
        plt.xlabel(r"$\eta$")
        plt.xlim(-8, 6)
        plt.ylim(-1.5, 2)
        plt.legend()
        plt.title(f"Solution F, λ = {self.lam}")
        plt.grid(True)
        sauver_figure(f"F_lambda_{self.lam}.png")

    def tracer_G(self):
        plt.figure()
        plt.plot(self.eta, self.G, label="G")
        plt.plot(self.eta, self.H, label="G'")
        plt.plot(self.eta, self.dH, label="G''")
        plt.xlabel(r"$\eta$")
        plt.xlim(-12, 6)
        plt.ylim(-1.5, 2)
        plt.legend()
        plt.title(f"Solution G, λ = {self.lam}, Pr = {self.Pr}")
        plt.grid(True)
        sauver_figure(f"G_lambda_{self.lam}_Pr_{self.Pr}.png")

    def tracer_champ_T(self, U1, nu1, L):
        """Champ T(x, y) : pour chaque x, y = delta(x) * (y/delta)(eta)."""
        x = np.linspace(L / 100, L, 200)
        delta = np.sqrt(2 * nu1 * x / U1)

        Y = np.concatenate(([-1000], self.Y, [1000]))
        T = np.concatenate(([self.T[0]], self.T, [self.T[-1]]))
        XX = np.outer(x, np.ones(len(Y)))
        YY = np.outer(delta, Y)
        TT = np.outer(np.ones(len(x)), T)

        plt.figure(figsize=(9, 4.5))
        cf = plt.contourf(XX * 100, YY * 1000, TT, levels=30, cmap="inferno")
        plt.colorbar(cf, label="T (K)")
        plt.ylim(-12 * delta[-1] * 1000, 5 * delta[-1] * 1000)
        plt.xlabel("x (cm)")
        plt.ylabel("y (mm)")
        plt.title(f"Champ de température, λ = {self.lam}"
                  f", U1 = {U1} m/s, T2 = {self.T2:.0f} K")
        sauver_figure(f"champ_T_lambda_{self.lam}.png")

    def sauver_profils(self):
        tableau = np.column_stack((self.eta, self.X[:, 0],
                                   self.X[:, 1], self.X[:, 2],
                                   self.G, self.H, self.dH, self.T, self.rho_s, self.Y))
        entete = "eta  F  F'  F''  G  G'  G''  T[K]  rho/rho1  y/delta"
        os.makedirs(DOSSIER, exist_ok=True)
        np.savetxt(f"{DOSSIER}/profils_lambda_{self.lam}_Pr_{self.Pr}.dat",
                   tableau,fmt="%12.6f", header=entete)


# ---------------------------------------------------------------
# Outils
# ---------------------------------------------------------------

def lire_parametres(nom="parametres.yaml"):
    """Lit le fichier YAML et renvoie un dictionnaire."""
    with open(nom, "r", encoding="utf-8") as fichier:
        p = yaml.safe_load(fichier)
    return p


def integrale(f, eta):
    """Primitive de f par trapèzes, qui vaut 0 en eta = 0."""
    I = np.zeros(len(eta))
    for i in range(len(eta) - 1):
        I[i + 1] = I[i] + 0.5 * (f[i] + f[i + 1]) * (eta[i + 1] - eta[i])
    return I - np.interp(0, eta, I)


def sauver_figure(nom):
    """Enregistre la figure courante dans Outputs
    (même dpi pour toutes) puis la ferme."""
    os.makedirs(DOSSIER, exist_ok=True)
    plt.tight_layout()
    plt.savefig(f"{DOSSIER}/{nom}", dpi=200)
    plt.close()


def afficher_tableau(lignes, entetes, titre, nom_fichier, floatfmt=".4f"):
    """Affiche le tableau à l'écran et l'enregistre (ascii) dans Outputs."""
    texte = tabulate(lignes, headers=entetes, floatfmt=floatfmt, tablefmt="simple_grid")
    print(f"\n{titre}\n{texte}")
    os.makedirs(DOSSIER, exist_ok=True)
    with open(f"{DOSSIER}/{nom_fichier}", "w", encoding="utf-8") as fichier:
        fichier.write(f"{titre}\n\n{texte}\n")


# ---------------------------------------------------------------
# Calcul pour tous les lambda
# ---------------------------------------------------------------

def etude_lambda(p, Pr, T2):

    n = p["numerique"]
    couches, r_eta, r_comp = [], [], []
    s = np.array(n["s0"])
    for lam in p["physique"]["lambdas"]:
        c = CoucheCisaillement(lam, Pr, p["physique"]["T1"], T2, n["L_haut"], n["L_zero"], n["h"])
        s = c.calculer(s, n["eps"], n["tol"], n["itmax"])
        couches.append(c)
        r_eta.append(c.ligne_eta)
        r_comp.append(c.ligne_comp)
    return couches, np.array(r_eta), np.array(r_comp)


# ---------------------------------------------------------------
# Options du menu
# ---------------------------------------------------------------

def sorties_F(couches, r_eta):
    """Option 1 : F, F', F'' et bornes choisies."""
    for c in couches:
        c.tracer_F()
    afficher_tableau(r_eta[:, [0, 1, 2, 3, 13, 14]],
                     ["λ", "F(0)", "F'(0)", "F''(0)", "η bas", "η haut"],
                     "Solution F en η = 0 et bornes du domaine", "F_eta_0.dat")

def sorties_post_F(r_eta):
    """Option 2 : grandeurs déduites de F (variable eta, densité constante)."""
    afficher_tableau(r_eta[:, [0, 4, 5, 6, 7, 8, 9]],
                     ["λ", "η_i", "u_i", "η_ω", "v(+∞)", "v(-∞)", "Cf (ρ = cte)"],
                     "Ligne de partage, épaisseur de vorticité, vitesses normales, frottement",
                     "post_F.dat")

def sorties_G(couches, r_eta):
    """Option 3 : G, G', G'' et Nu."""
    for c in couches:
        c.tracer_G()
    afficher_tableau(r_eta[:, [0, 10, 11, 12]], ["λ", "G'(0)", "G''(0)", "Nu (ρ = cte)"],
                     "Solution G en η = 0 et nombre de Nusselt", "G_eta_0.dat")

def sorties_compressible(couches, r_comp, p):
    """Option 4 : profils u/U1 et T/T1 en y/delta (tous les lambda), champ T(x, y)."""
    for nom, etiquette in [("u", r"$u/U_1$"), ("T", r"$T/T_1$")]:
        plt.figure()
        for c in couches:
            valeurs = c.X[:, 1] if nom == "u" else c.T / c.T1
            plt.plot(valeurs, c.Y, label=f"λ = {c.lam}")
        plt.xlabel(etiquette)
        plt.ylabel(r"$y/\delta$")
        plt.ylim(-20, 8)
        plt.legend()
        plt.grid(True)
        plt.title(f"Profils compressibles, T1 = {couches[0].T1:.0f} K, T2 = {couches[0].T2:.0f} K")
        sauver_figure(f"{nom}_compressible.png")
    d = p["dimensionnel"]
    for c in couches:
        c.tracer_champ_T(d["U1"], d["nu1"], d["L"])
        c.sauver_profils()
    afficher_tableau(r_comp, ["λ", "ρ*(0)", "δω/δ", "Cf", "Nu"],
                     "Grandeurs compressibles", "compressible.dat")

def etude_T2(p, liste_T2):
    """Option 5 : grandeurs compressibles pour plusieurs T2."""
    lignes = []
    for T2 in liste_T2:
        couches, r_eta, r_comp = etude_lambda(p, p["physique"]["Pr"], T2)
        for ligne in r_comp:
            lignes.append([T2, *ligne])
    afficher_tableau(lignes, ["T2 (K)", "λ", "ρ*(0)", "δω/δ", "Cf", "Nu"],
                     "Étude paramétrique en T2", "etude_T2.dat")


# ---------------------------------------------------------------
# Validation
# ---------------------------------------------------------------

def lire_valeurs_prof(nom):
    valeurs = {}
    with open(nom, "r", encoding="utf-8") as fichier:
        for ligne in fichier:
            if "||" not in ligne:
                continue
            nom_g, reste = ligne.split("||")
            nombres = [float(v) for v in reste.split("|") if v.strip() != ""]
            valeurs[nom_g.strip()] = nombres
    return valeurs


def comparer_prof(resultats, nom_fichier_prof):
    """Nos valeurs converties dans la convention
    du prof (eta_prof = sqrt(2) eta)."""
    prof = lire_valeurs_prof(nom_fichier_prof)
    r2 = np.sqrt(2)
    correspondance = [("f0", 1, r2), ("u0", 2, 1), ("du0", 3, 1 / r2),
                      ("etai", 4, r2), ("ui", 5, 1), ("etaOmega", 6, r2),
                      ("v{+inf}", 7, 1), ("v{-inf}", 8, 1)]
    lignes = []
    for i in range(len(resultats)):
        for nom_g, col, facteur in correspondance:
            nous = resultats[i, col] * facteur
            ref = prof[nom_g][i]
            lignes.append([resultats[i, 0], nom_g, nous, ref, nous - ref])
    afficher_tableau(lignes, ["λ", "grandeur", "code", "prof", "écart"],
                     "Comparaison avec les valeurs du prof",
                     "validation_prof.dat")


def comparer_schlichting(resultats):
    """Schlichting & Gersten, tableaux 7.1 et 22.2 :
    lambda = 0, Pr = 0,5 (comparaison directe)."""
    ligne = resultats[0]
    livre = [("F(0)", 1, 0.2392), ("F'(0)", 2, 0.6914), ("F''(0)", 3, 0.2704),
             ("eta_i", 4, -0.3740), ("u_i", 5, 0.5872),
             ("G'(0)", 10, 0.6937), ("G''(0)", 11, 0.1823)]
    lignes = []
    for nom, col, ref in livre:
        lignes.append([nom, ligne[col], ref, ligne[col] - ref])
    afficher_tableau(lignes, ["grandeur", "code", "livre", "écart"],
                     "Comparaison avec Schlichting (λ = 0)",
                     "validation_schlichting.dat")


def verifier_reynolds(resultats):
    """À Pr = 1 : G'(0) = (F'(0) - lambda)/(1 - lambda) et Nu = -Cf/(1 - lambda)."""
    lignes = []
    for ligne in resultats:
        lam, u0, Cf, Gp0, Nu = ligne[0], ligne[2], ligne[9], ligne[10], ligne[12]
        lignes.append([lam, Gp0, (u0 - lam) / (1 - lam), Nu, -Cf / (1 - lam)])
    afficher_tableau(lignes, ["λ", "G'(0)", "attendu", "Nu", "attendu"],
                     "Analogie de Reynolds (Pr = 1)",
                     "validation_reynolds.dat")
