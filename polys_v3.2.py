#Projet MIDL Polynômes


# -*- coding: utf-8 -*-

import cmath
import copy
from math import *
import numpy
from fractions import Fraction

# =============================================================================
# Représentation des données
# =============================================================================

"""
polynome :
    liste de coefficients
    /!\ le terme de plus haut degré est le plus à droite
    [1,0,56,1,0,0]  <=>  1 + 56.X**2 + X**3
"""

#Polynômes réels
p1 = [3,0,9,0,0,1,2]
p2 = [0,2,3,5,0,0,0]
p3 = [0,0,1]
p4 = [30,0,4,8,0,0,0,0]
p5 = [56,4,5,3]
p6 = [1,3,3,1]
p7 = [1,1,1,1]
p8 = [0,40000,0,1]
p_einstein = [-12,-14,0,1]

#Polynômes complexes
pc1 = [2,complex(1,2),complex(5,6)]
pc2 = [complex(3,1),complex(-6,-7),complex(0,3),complex(2,0)]

# =============================================================================
#  TESTS DES FONCTIONS
# =============================================================================

def test_add():
    assert(add(0,1) == 1)
    assert(add(1,1) == 2)

def test_sub():
    assert(sub(0,1) == -1)
    assert(sub(1,1) == 0)


def test_deg():
    assert(deg([]) == 0)
    assert(deg([0,0,0]) == 0), "le polynône nul peut avoir plusieurs 0"
    assert(deg([1,0,1,1,0,0]) == 3), "des 0 devant n'affectent pas le degré"
    assert(deg([1,0,1,1]) == 3)

def test_coeff():
    assert(coeff([]) == 0)
    assert(coeff([0,0]) == 0)
    assert(coeff([0,1,1,0]) == 1), "0 ne peut pas être un coefficient dominant"
    assert(coeff([0,1,1]) == 1), "les termes de plus haut degré sont à la fin de la liste"

def test_is_null_poly():
    assert(is_null_poly([]))
    assert(is_null_poly([0,0])), "le polynome nul peut avoir plusieurs 0"
    assert(not(is_null_poly([0,1,0])))

def test_add_poly():
    assert(add_poly([0,0,1,1], [1]) == [1,0,1,1])
    assert(add_poly([0,1,1], [1,0,0,1]) == [1,1,1,1])
    assert(add_poly([0,1,1], [1,0,1,1]) == [1,1,2,1])
    assert(add_poly([0,1,1,0], [1,0,1,1,0]) == [1,1,2,1,0])

def test_mult_par_monome():
    assert(mult_par_monome([0,0,1,0,0],1,[0,1,0,0,1,1,0,0]) == [0,0,0,1,0,0,1,1,0,0])
    assert(mult_par_monome([0,0,0,0,1],3,[0,1,1,1]) == [0,0,0,0,0,3,3,3]), "multiplier par le coefficient chaque terme"

def test_mult_poly():
    assert(mult_poly([0,1,1], [1,0,1,1]) == [0,1,1,1,2,1])

def test_puissance_poly():
    poly = [3,0,9,0,0,1,2]
    assert(puissance_poly(poly,2) == [9,0,54,0,81,6,12,18,36,0,1,4,4])
    assert(puissance_poly([1,1],3) == [1,3,3,1])

def test_eucl_poly():
    assert((eucl_poly([0,1,1,0,1],[1,1])) == ([-1,2,-1,1], [1,0,0,0,0]))
    assert(div_poly([0,1,1,0,1], [1,1]) == [-1,2,-1,1])
    assert(mod_poly([0,1,1,0,1], [1,1]) == [1,0,0,0,0])

def test_derivee(poly_reel,poly_complexe):
    assert(derivee([9]) == [0]), "dérivée d'un constante nulle"
    assert(derivee(poly_reel) == [2,6,15])
    assert(derivee(poly_complexe) == [complex(-6,-7),complex(0,6),complex(6,0)])

def test_calcul(poly_reel,poly_comp):
    assert(calcul(poly_reel,1) == 10)
    assert(calcul(poly_reel,0) == 0),"si X=0, égal au premier terme (sans X)"
    assert(calcul(poly_comp,complex(1,2)) == complex(-23,-31))
    assert(calcul(poly_comp,0) == complex(3,1)),"si X=0, égal au premier terme (sans X)"


def test_delta_2():
    pass

def test_racines():
    pass

def test_racines_cubiques(p_einstein):
    """attention int(...) n'arrondit pas mais donne la partie entière"""
    racine1_cube = puissance_poly([1,1],3)
    r1,r2,r3 = racines_cubiques(racine1_cube)
    ##print(int(r1),r1,int(r2),r2,int(r3),r3)
    assert({r1,r2,r3} == {-1,-1,-1})
    #Cas où il n'y a pas de terme de degré 2
    r1,r2,r3 = racines_cubiques(p_einstein)
    assert({int(r1),int(r2),int(r3)} == {-3,0,4})
    r1,r2,r3 = racines_cubiques([-100,-96,1])
    assert({int(r1),int(r2),int(r3)} == {-9,-1,9})
    #Cas où il y a des termes de degré 2
    r1,r2,r3 = racines_cubiques([-256,-64,8,1])
    assert({int(r1),int(r2),int(r3)} == {-11,-3,6}) #on prend un ensemble car on ne sait pas dans quel ordre les racines apparaissent

def approx_eq(n1,n2):
    #Approximation grossière. Bien meilleure que la partie entière cependant
    if abs(n1-n2) < 10**-1:
        return True
    return False
def approx_eq_liste(l1,l2):
    #Regarde si deux listes de nombres sont égales à 10**-2 près (une fois mises en ordre croissant)
    if(len(l1) != len(l2)): return False
    l1.sort()
    l2.sort()
    
    for i in range(len(l1)):
        if not approx_eq(l1[i],l2[i]):
            #print(l1[i],l2[i])
            return False
    return True

def test_racines_quartiques():
    p1 = [-6,-13,-7,1,1]
    r1,r2,r3,r4 = racines_quartiques(p1)
    #print(r1,r2,r3,r4)
    assert(approx_eq_liste([r1.real,r2.real,r3.real,r4.real],[-1,-1,3,-2]))
    p2 = [-10000,0,0,0,1]
    r1,r2,r3,r4 = racines_quartiques(p2)
    #print("Les racines sont",r1,r2,r3,r4)
    assert(approx_eq_liste([r1.real,r2.real,r3.real,r4.real],[0,-10,10,0]))

    p3 = [0,1,1,1,1]
    r1,r2,r3,r4 = racines_quartiques(p3)
    #print("Les racines sont",r1,r2,r3,r4)
    assert(approx_eq_liste([r1.real,r2.real,r3.real,r4.real],[-1,0,0,0]))
    assert(approx_eq_liste([r1.imag,r2.imag,r3.imag,r4.imag],[1,-1,0,0]))

def tests():
    """effectue l'ensemble des tests"""
    poly_reel = [0,2,3,5,0,0,0]
    poly_complexe = [complex(3,1),complex(-6,-7),complex(0,3),complex(2,0)]
    p_einstein = [-12,-14,0,1]
    test_deg()
    test_coeff()
    test_is_null_poly()
    print("Tests des fonctions d'accès aux polynômes : OK")
    test_add()
    test_sub()
    test_add_poly()
    print("Tests relatifs à l'addition / soustraction : OK")
    test_mult_par_monome()
    test_mult_poly()
    test_eucl_poly()
    test_puissance_poly()
    print("Tests relatifs à la multiplication / division : OK")
    test_derivee(poly_reel,poly_complexe)
    test_calcul(poly_reel,poly_complexe)
    print("Autres Tests : OK")
    test_delta_2()
    test_racines()
    print("Tests des polynômes de degré 2 : OK")
    test_racines_cubiques(p_einstein)
    print("Tests des polynômes de degré 3 : OK")
    test_racines_quartiques()
    print("Tests des polynômes de degré 4 : OK")



# =============================================================================
#  DEFINITION DES FONCTIONS
# =============================================================================

# =============================================================================
# Représentation graphique
# =============================================================================

def affiche_poly(poly):
    degre = 0
    premier_non_nul = True
    for coefficient in poly:
        if coefficient == 0:
            print("      ",end=" ")
        else:
            if (not(premier_non_nul)):
                print("+ ",end="")
            premier_non_nul = False
            if degre == 0:
                print(coefficient,"     ",end=" ")
            elif degre == 1:
                print(str(coefficient)+".X   ",end=" ")
            else:
                if coefficient == 1: print("  X**"+str(degre),end=" ") ##str supprime l'espace en trop
                else: print(str(coefficient)+".X**"+str(degre),end=" ")
        degre+=1

def afficher_calculs(poly1,poly2,symbole):
    """affiche les 2 polynômes, le symbole de l'opération et le résultat.
    /!\ pour l'instant cette fonction n'affiche que l'addition
    """
    print("  ",end="")
    affiche_poly(poly1)
    print("")
    print(symbole+" ",end="")
    affiche_poly(poly2)
    print("")
    print(" _"+"_"*9*(max(deg(poly1),deg(poly2))+1))
    result = add_poly(poly1,poly2)
    print("  ",end="")
    affiche_poly(result)
    return result


# =============================================================================
#  Fonctions Unitaires
# =============================================================================

#Addition

def add(n0,n1):
    return n0+n1

def sub(n0,n1):
    return n0-n1

#Degré d'un polynome

def deg(poly):
    degre=len(poly)-1
    while degre>=0 and poly[degre]==0:
        degre-=1
    if degre==-1:
        return 0
    return degre


#Coefficient dominant du polynôme (le premier différent de 0)

def coeff(poly):
    l=len(poly)-1
    while l>0:
        if poly[l]:
            return 1
        l-=1
    return 0

#p est-il le polynôme nul ?

def is_null_poly(p):
    if coeff(p)==0 and deg(p)==0:
        return True
    return False


#Addition de deux polynômes

def add_poly(p1,p2):
    deg1,deg2=deg(p1),deg(p2)
    if p1==[]:
        return(p2)
    if p2==[]:
        return(p1)
    p_somme=[]
    nb_coeffs_max=len(p1)
    if nb_coeffs_max<len(p2):
        nb_coeffs_max=len(p2)
    for i_coeff in range(nb_coeffs_max):
        if i_coeff<=deg1 and i_coeff<=deg2:
            p_somme.append(add(p1[i_coeff],p2[i_coeff]))
        elif i_coeff<=deg1:
            p_somme.append(p1[i_coeff])
        elif i_coeff<=deg2:
            p_somme.append(p2[i_coeff])
        else:
            p_somme.append(0)
    return(p_somme)

#Soustraction de 2 polynômes

def diff_poly(p1,p2):
    deg1,deg2=deg(p1),deg(p2)
    if p2==[]:
        return(p1)
    p_somme=[]
    nb_coeffs_max=len(p1)
    if nb_coeffs_max<len(p2):
        nb_coeffs_max=len(p2)
    for i_coeff in range(nb_coeffs_max):
        if i_coeff<=deg1 and i_coeff<=deg2:
            p_somme.append(sub(p1[i_coeff],p2[i_coeff]))
        elif i_coeff<=deg1:
            p_somme.append(p1[i_coeff])
        elif i_coeff<=deg2:
            p_somme.append(-p2[i_coeff])
        else:
            p_somme.append(0)
    return(p_somme)


#Multiplier 2 polynomes

def mult_par_monome(monome,coeff,p):
    poly=p[:]
    degre=deg(monome)
    while degre>0:
        poly.insert(0,0)
        degre-=1
    for i_coeff in range(len(poly)):
        poly[i_coeff]=poly[i_coeff]*coeff
    return poly

def mult_poly(p1,p2):
    poly_inter=[] #polynome intermédiaire pour les calculs
    poly_final=[]
    monome=[] #liste de 0
    for coeff in p1:
        if coeff:
            poly_inter=mult_par_monome(monome+[1],coeff,p2)
            poly_final=add_poly(poly_inter,poly_final)
        monome.append(0)
    return poly_final

def puissance_poly(p,n):
    result = [1]
    for i in range(n):
        ##print(result,p)
        result = mult_poly(result,p)
        ##print("multiplication en cours",i,result)
    return result



#Division euclidienne de 2 polynomes

def creer_monome(degre):
    poly=[]
    for i in range(degre):
        poly.append(0)
    poly.append(1)
    return(poly)

def eucl_poly(p1,p2):
    #Cassé en dehors de Z/2Z. Tu peux supprimer si tu veux la partie en """"""
    """
    q=[]
    r=p1
    d=deg(p2)
    c=coeff(p2)
    print(p1,p2)
    while not(is_null_poly(r)) and deg(r)>=d:
        s=mult_par_monome(creer_monome(deg(r)-d),int(coeff(r)/c),[1])
        q=add_poly(q,s)
        r=diff_poly(r,mult_poly(s,p2))
        print(s,q,r)
    return (q,r)
    """
    q,r = numpy.polydiv(numpy.array(p1[::-1]), numpy.array(p2[::-1]))
    r = list(r[::-1])
    while(len(r) < max(len(p1),len(p2))):
        r.append(0)
    return list(q[::-1]),r


#Quotient et reste de la division euclidienne de 2 polynômes

def div_poly(p1,p2):
    (q,r)=eucl_poly(p1,p2)
    return q

def mod_poly(p1,p2):
    (q,r)=eucl_poly(p1,p2)
    return r

#Dérive le polynôme

def derivee(poly):
    """Dérive poly"""
    poly_d = []
    for i in range(deg(poly)+1):
        poly_d.append(poly[i]*i)
    if (len(poly) > 1):
        poly_d.pop(0)
    return poly_d


def calcul(poly,valeur):
    """remplace X par valeur dans le polynôme"""
    result = 0
    for i in range(deg(poly)+1):
        result+=poly[i]*(valeur**i)
    return result


# =============================================================================
#  Equations quadratiques
# =============================================================================

def delta_2(p):
    """discriminant d'un polynôme de degré 2"""
    a = p[2]
    b = p[1]
    c = p[0]
    return b*b-4*a*c

def racines(poly):
    """renvoie les racines de poly"""
    discriminant=delta_2(poly)
    a=poly[2]
    b=poly[1]
    x1=(-b-cmath.sqrt(discriminant))/(2*a)
    x2=(-b+cmath.sqrt(discriminant))/(2*a)
    if (x1.imag == 0 and x2.imag == 0):
        return (x1.real,x2.real)
    else:
        return (x1,x2)



# =============================================================================
#  Equations cubiques
# =============================================================================


def racines_cubiques(p):
    """Méthode de Cardan améliorée avec les nombres complexes"""
    #print("Polynôme :")
    #affiche_poly(p)
    #print("")
    coeff_dominant = coeff(p)
    if coeff_dominant != 1:
        for i in range(len(p)):
            p[i] = p[i]/coeff_dominant #le coefficient dominant devient 1
    if p[2] != 0:
        transfo_Tchirnhaus = [(-1/3)*p[2],1]
        p_deg3 = puissance_poly(transfo_Tchirnhaus,3)
        p_deg2 = puissance_poly(transfo_Tchirnhaus,2)
        #Suppression des termes de degré 2
        p_result = add_poly(add_poly(p_deg3,mult_poly(p_deg2,[p[2]])),add_poly(mult_poly(transfo_Tchirnhaus,[p[1]]),[p[0]]))
    else:
        p_result = p
    #print("p result :")
    #affiche_poly(p_result)
    #print("")

    #Trouver une racine
    p1 = p_result[1]
    p2 = p_result[0]
    ##print( 4*(p1**3)+27*(p2**2), p1,p2,4*(p1**3),27*(p2**2))
    if 4*(p1**3)+27*(p2**2) >= 0:
        ##print("if")
        r_delta = sqrt((4*(p1**3)+27*(p2**2))/27)
        ##print("rrrrr",r_delta)
        a = (-p2-r_delta)/2
        b = (-p2+r_delta)/2
    else:
        ##print("else")
        r_delta = cmath.sqrt((4*(p1**3)+27*(p2**2))/27)
        a = (-p2-r_delta)/2
        b = (-p2+r_delta)/2
    ##print("A ET B",a,b)
    ##print(-1/((-a)**3))
    ##print("a :",a)
    ##print("type de (a**Fraction(1,3)) :",type(a**Fraction(1,3)))
    ##print(a**Fraction(1,3), b**Fraction(1,3))
    if a.imag == 0:
        ##print("reeeeel")
        x1 = -(-a)**Fraction(1,3) + b**Fraction(1,3)
    else:
        x1 = a**Fraction(1,3) + b**Fraction(1,3)
    if (x1.imag == 0):
        x1 = x1.real
    ##print(x1,type(x1))

    #Déterminer les autres
    poly_reste = div_poly(p_result,[-x1,1])
    result = calcul(p_result,x1)
    #affiche_poly(p_result)
    #print(x1)
    assert(complex(int(result.real),int(result.imag)) == complex(0,0)) #x1 est une racine de p_result
    #int permet d'arrondir
    x2,x3 = racines(poly_reste)
    return (x1-p[2]/3,x2-p[2]/3,x3-p[2]/3)



# =============================================================================
#  Equations quartiques
# =============================================================================

def racines_quartiques(p):
    #Conversion en un polynôme de forme y^4 + py^2 + qy + r
    coeff_dominant = coeff(p)
    for i in range(len(p)):
        p[i] = p[i]/coeff_dominant #le coefficient dominant devient 1
    #x = y-p[-2]/4
    transfo_depressed = [(-1/4)*p[-2],1]
    p_deg4 = puissance_poly(transfo_depressed,4)
    p_deg3 = puissance_poly(transfo_depressed,3)
    p_deg2 = puissance_poly(transfo_depressed,2)
    p_result = add_poly(mult_poly(p_deg4,[p[4]]),add_poly(add_poly(mult_poly(p_deg3,[p[3]]),mult_poly(p_deg2,[p[2]])),add_poly(mult_poly(transfo_depressed,[p[1]]),[p[0]])))
    """
    On a donc y^4 + py^2 + qy + r
    Et on cherche (y^2 + ay + b)(y^2 + cy + d)
    c = -a
    b = a^2/2 + p/2 - q/2a
    d = a^2/2 + p/2 + q/2a
    On cherche donc EQ1 : a^6 + 2pa^4 + (p^2-4r)a^2 - q^2 = 0
    Soit b = a^2
    On a alors EQ2 : b^3 + 2pb^2 + (p^2-4r)b - q^2 = 0
    """
    #print("EQ2 : ",[-p_result[1]**2,p_result[-3]**2-4*p_result[0],2*p_result[-3],1])
    EQ2 = racines_cubiques([-p_result[1]**2,p_result[-3]**2-4*p_result[0],2*p_result[-3],1])
    a = cmath.sqrt(EQ2[1])
    #On a donc
    #print(a)
    P1 = [a**2/2 + p_result[-3]/2 - (p_result[1]/(2*a)),a,1]
    P2 = [a**2/2 + p_result[-3]/2 + (p_result[1]/(2*a)),-a,1]
    #Avec P1*P2 = p
    x1,x2 = racines(P1)
    x3,x4 = racines(P2)
    #print(P1,P2)
    x1 = x1 - p[3]/4
    x2 = x2 - p[3]/4
    x3 = x3 - p[3]/4
    x4 = x4 - p[3]/4
    return(x1,x2,x3,x4)



# =============================================================================
#  Navigation
# =============================================================================

def menu_principal():
    print("Menu :")
    print("0 : Exécuter l'ensemble des tests")
    print("1 : A voir")
    print("2 : A voir")
    print("3 : A voir")
    print("4 : A voir")
    print("")

def choix_exec(msg,possibilites):
    """renvoie un choix d'exécution valide"""
    partie=input(msg)
    while partie not in possibilites:
        partie=input(msg)
    return partie


# =============================================================================
#  Exécution des tests
# =============================================================================

def exec_test():
    print(" ----------------- Lancement des tests ----------------- ")
    print("")
    tests()
    print("")
    print("--------- Tous les tests ont bien été exécutés --------- ")


# =============================================================================
#  Programme principal
# =============================================================================

menu_principal()

choix = choix_exec("Votre choix (0, 1, 2, 3 ou 4) : ",('0','1','2','3','4'))

if choix == '0': #Exécution de l'ensemble des tests
    exec_test()
else:
    r1,r2,r3 = racines_cubiques(p6)
    print(r1,r2,r3)
    r1,r2,r3 = racines_cubiques(p_einstein)
    print(r1,r2,r3)
    r1,r2,r3 = racines_cubiques(p8)
    print(r1,r2,r3)
    r1,r2,r3 = racines_cubiques(p7)
    print(r1,r2,r3)
    racines_quartiques([20,-16,-9,4,1])