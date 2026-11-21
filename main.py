import couche as cc

p = cc.lire_parametres()
Pr = p["physique"]["Pr"]
T2 = p["physique"]["T2"]

print("1 : Résoudre et tracer F (F, F', F'')")
print("2 : Post-traitement de F (ligne de partage, δω, vitesses normales, Cf)")
print("3 : Résoudre et tracer G (G, G', G'') et Nu")
print("4 : Retour au compressible (profils en y, champ T(x,y), δω, Cf, Nu)")
print("5 : Étude paramétrique en T2")
print("6 : Validation : valeurs du prof et Schlichting")
print("7 : Validation : analogie de Reynolds (Pr = 1)")
print("8 : Tout (options 1 à 4)")
choix = input("Choisir une option : ")

if choix in ["1", "2", "3", "4", "6", "8"]:
    couches, r_eta, r_comp = cc.etude_lambda(p, Pr, T2)

if choix in ["1", "8"]:
    cc.sorties_F(couches, r_eta)

if choix in ["2", "8"]:
    cc.sorties_post_F(r_eta)

if choix in ["3", "8"]:
    cc.sorties_G(couches, r_eta)

if choix in ["4", "8"]:
    cc.sorties_compressible(couches, r_comp, p)

if choix in ["5", "8"]:
    cc.etude_T2(p, [300, 450, 600, 900])

if choix in ["6", "8"]:
    cc.comparer_prof(r_eta, "Valeurs.dat.txt")
    cc.comparer_schlichting(r_eta)

if choix in ["7", "8"]:
    couches, r_eta, r_comp = cc.etude_lambda(p, 1.0, T2)
    cc.verifier_reynolds(r_eta)

if choix not in ["1", "2", "3", "4", "5", "6", "7", "8"]:
    print("Option inconnue")
