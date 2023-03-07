#Projet de Mathématiques _ Implémenter l'algorithme de Waring
#Ducros Alexandre, Bourkdeche Melissa, Cauquil Diane L2 MIDL

# -*- coding: utf-8 -*-

## /!\ IMPORTANT => On part du principe que le polynôme entré est symétrique

import copy
import sys
sys.setrecursionlimit(1000000)
# =============================================================================
# Règles de codage
# =============================================================================

""" Tous les noms des variables et des fonctions sont écrits en français.
    Les noms des fonctions sont écrits en Snake... (avec des underscores / sans majuscules)
    Aucune fonction ne dépasse les 50 lignes.
    Chaque fonction est nécessaire et contient de la docstring.
    Les nombres (coefficients des polynômes) sont des float."""


# =============================================================================
# Représentation des données (classes)
# =============================================================================

class Objet:
    pass

class Nombre(Objet):
    def __init__(self,value):
        self.value = value
    def to_str(self,debug=False):
        return str(self.value)

class Inconnue(Objet):
    def __init__(self, name):
        self.name = name
        self.quantity = 1.0
        self.puissance = 1.0

    def to_str(self,debug=False):
        s = ""
        if self.quantity != 1.0:
            s += str(self.quantity) + "*"
        s += self.name
        if self.puissance != 1.0:
            s += "^"+str(self.puissance)
        return s

class Op(Objet):
    def __init__(self,*elems):
        self.elems = list(elems)
        self.op = ""
    def to_str(self, debug = False):
        if debug:
            s = str(type(self).__name__) + "("
            for x in range(len(self.elems)):
                s += self.elems[x].to_str(debug=True)
                if x != len(self.elems)-1:
                    s += ","
            s += ")"
            return s
        s = "("
        for x in range(len(self.elems)):
            s += self.elems[x].to_str()
            if x != len(self.elems)-1:
                s += self.op
        s+=")"
        return s

class Add(Op):
    def __init__(self,*elems):
        Op.__init__(self,*elems)
        self.op = "+"
    def resultat(self):
        somme = 0
        for elem in self.elems:
            somme+=elem.value
        return somme

class Mult(Op):
    def __init__(self,*elems):
        Op.__init__(self,*elems)
        self.op = "*"
    def resultat(self):
        produit = 1
        for elem in self.elems:
            produit *= elem.value
        return produit

class Div(Op):
    def __init__(self,left,right):
        Op.__init__(self,right,left)
        self.op = "/"
    def resultat(self):
        return self.elems[0].value / self.elems[1].value

class Sub(Op):
    def __init__(self,*elems):
        Op.__init__(self,*elems)
        self.op = "-"
        if len(self.elems) > 0:
            self.elems = self.elems[::-1]
    def resultat(self):
        orig = self.elems[0].value
        for x in range(1,len(self.elems)):
            orig -= self.elems[x].value
        return orig

class Puis(Op):
    def __init__(self, left,right):
        Op.__init__(self, right,left)
        self.op = "^"
    def resultat(self):
        return self.elems[0]**self.elems[1]


# =============================================================================
#  DEFINITION DES FONCTIONS
# =============================================================================

# =============================================================================
# Fonctions Unitaires pour trouver des informations (lecture)
# =============================================================================

def precedence(token):
    """renvoie la priorité de l'opérateur token (un symbole de type str)
    token peut être une chaîne de plusieurs caractères (nom d'une inconnue)"""
    precedence_table = {
        "(" : 13,
        ")" : 13,
        "^" : 13,
        "*" : 10,
        "/" : 10,
        "+" : 9,
        "-" : 9
    }
    return precedence_table[token]

def left_associative(token):
    """renvoie true si l'opérateur token est associatif à gauche"""
    left_associative_table = {
        "+" : True,
        "*" : True,
        "/" : True,
        "-" : True,
        "^" : False
    }
    return left_associative_table[token]

def est_nombres(liste):
    """Renvoie True si liste est une liste constituée uniquement d'éléments de type Nombre"""
    for elem in liste:
        if type(elem) != Nombre:
            return False
    return True

def est_expression(chaine):
    li_toks = get_tokens(chaine)
    for elem in li_toks:
        if type(elem) == str and elem not in ("+","*","/","-","^","(",")"):
            return False
    return True

#Monômes et coefficients
def degre_monome(monome):
    """renvoie le degré du monome (un arbre), les inconnues
    sont séparées par des * ou des ^
    /!\ les ^ sont encadrés des 2 côtés par des parenthèses"""
    m = str_to_tree(monome)
    deg = 0
    if type(m) == Nombre:
        return 0
    if type(m) == Inconnue:
        return 1
    for terme in m.elems:
        deg+=terme.puissance
    return deg

##print(degre_monome("x1^2*x2"))
##print(degre_monome("x3"))
##print(degre_monome("1"))

def terme_dominant(poly):
    #print("on appelle terme dominant")
    inconnues = str_to_dico(poly)
    deg_max = 0
    monome = '0'
    for inconnue in inconnues:
        deg = degre_monome(inconnue)
        if deg_max < deg:
            deg_max = deg
            monome = inconnue
    return monome

###print("Dominus rex",terme_dominant("x1^2*x2+x3"))

def coeff_de(inconnue, poly):
    #print("on appelle coeff de")
    """renvoie le coefficient devant l'inconnue dans poly (0 si le terme n'apparaît pas"""
    dico = str_to_dico(poly)
    if inconnue in dico:
        return float(dico[inconnue])
    return 0.0

def list_pow(monome,nb_inconnue):
    #print("on appelle list pow")
    """renvoie la liste des puissances du monôme (un arbre) dans l'ordre
    d'apparition des inconnues"""
    li_puis = []
    if type(monome) == Inconnue:
        li_puis=[monome.puissance]
    else:
        for inconnue in monome.elems:
            if type(inconnue) == Inconnue:
                li_puis.append(inconnue.puissance)

    while len(li_puis) <= nb_inconnue: #il faut un nombre pair de numéros
        li_puis.append(0.0)

    return li_puis

def noms_inconnues(poly):
    #print("on appelle noms inconnues")
    """renvoie la liste des noms des inconnues"""
    arbre = str_to_tree(poly)

    #c'est un polynôme donc des add de mult ou d'inconnues
    inconnues = set()
    for elem in arbre.elems:
        if type(elem) == Inconnue:
            inconnues.add(elem.name)
            continue
        if type(elem) == Nombre: continue

        for sub_elem in elem.elems:
            if type(sub_elem) == Inconnue:
                inconnues.add(sub_elem.name)
    return list(inconnues)

#Polynômes symétriques élémentaires
def parties(elem_list):
    #print("on appelle parties")
    liste = set()
    def parties_de(elems,cur):
        liste.add(tuple(cur))
        if len(elems) == 1:
            liste.add(tuple(cur+[elems[0]]))
            return
        parties_de(elems[:-1],cur + [elems[-1]]) #on prend
        parties_de(elems[:-1],cur) #on ne prend pas
    parties_de(elem_list,[])
    return list(liste)

def sigma(nb_inconnues, n):
    #print("on appelle sigma")
    """Renvoie une chaine de caractères représentant le n-ième polynôme symétrique
    élémentaire"""
    poly = ""
    inconnues = []
    for x in range(nb_inconnues):
        inconnues.append("x" + str(x+1))
    all_parties = parties(inconnues)
    liste = []
    for elem in all_parties:
        if len(elem) == n:
            liste.append(elem)
    for y in range(len(liste)):
        elem = liste[y]
        sub_poly = ""
        for x in range(len(elem)):
            sub_poly += elem[x]
            if x != len(elem)-1:
                sub_poly += "*"
        if y != len(liste)-1:
            sub_poly += " + "
        poly += sub_poly
    return poly

def jolis_sigmas(poly):
    sig = "\u03C3"
    return poly.replace("sigma",sig)

def remplacer_par(poly1, inconnue1, poly2):
    #print("on appelle remplacer par")
    """remplace inconnue1 par le poly1 dans le poly2 (tous sont des str)"""
    poly3 = poly2.replace(inconnue1,'('+poly1+')')
    return poly3

def cote_a_cote(c1,c2,c3,liste):
    #print("on appelle cote a cote")
    """renvoie True si dans liste on a : c2,c1 ou c1,c3"""
    i=0
    while (i+1<len(liste)):
        if ((liste[i]==c1 and liste[i-1]==c2) or (liste[i]==c1 and liste[i+1]==c3)):
            return True
        i+=1
    return False

#liste=['3.0','aaa',3.0,'*', 1, '*',')','*',5]
##print("COTE A COTE",cote_a_cote('*',')','(',liste))


# =============================================================================
# Fonctions Unitaires pour changer de représentation
# =============================================================================

def add_dans_dico(inconnue,coeff,dico):
    if inconnue in dico:
        dico[inconnue]+=coeff
    else:
        dico[inconnue]=coeff


# =============================================================================
# Fonctions Unitaires pour effectuer des calculs
# =============================================================================

#Soustraction
def fusion_sub(tokens):
    #print("on appelle fusion sub")
    """Fusionne les symboles - qui sont seuls avec les coefficients.
    Entrée / Sortie : une liste de tokens (les chiffres sont des float)"""
    i_tok = 0
    operateurs = ['(','+','*','^','/','-'] #')' n'en fait pas partie puisque (3)-2 a 2 opérandes
    while (i_tok<len(tokens)):
        if (tokens[i_tok] == '-' and (i_tok == 0 or (tokens[i_tok-1] in operateurs))) and type(tokens[i_tok+1]) == float: #c'est un moins solitaire
            tokens[i_tok] = -tokens[i_tok+1] #on le fusionne avec le chiffre qui suit
            del tokens[i_tok+1]
        elif tokens[i_tok] == '-' and type(tokens[i_tok+1]) == str and tokens[i_tok+1] not in ['+','*','^','/','-',')']: #j'ai recopié la liste pour rajouter le ) en +
            #on a donc 3-x1
            #on veut avoir 3 + -1*x1
            tokens[i_tok] = '+'
            tokens.insert(i_tok+1,-1.0)
            tokens.insert(i_tok+2,"*")
        elif tokens[i_tok] == '-' and type(tokens[i_tok+1]) == float: #c'est un moins encore plus solitaire (3-2-1 par exemple)
            tokens[i_tok] = '+'
            tokens[i_tok+1] *= -1.0
        i_tok+=1
    return tokens

#Addition
def compresser_add(arbre):
    i = 0
    while(i < len(arbre.elems)):
        j = 0
        trouve = False
        while(j < i and not trouve):
            if type(arbre.elems[j]) == type(arbre.elems[i]) == Nombre:
                arbre.elems[j].value += arbre.elems[i].value
                del arbre.elems[i]
                trouve = True
            j += 1
        if not trouve:
            i += 1
    return arbre

#Multiplication
def compresser_mult(arbre):
    for _ in range(10):
        i = 0
        while(i < len(arbre.elems)):
            j = 0
            trouve = False
            while(j < i and not trouve):
                if type(arbre.elems[j]) == Inconnue and type(arbre.elems[i]) == Nombre: #pour enlever le double if on fait ça
                    arbre.elems[j],arbre.elems[i] = arbre.elems[i],arbre.elems[j]

                if type(arbre.elems[j]) == type(arbre.elems[i]):
                    trouve = True
                    if type(arbre.elems[j]) == Nombre:
                        arbre.elems[j].value *= arbre.elems[i].value
                        del arbre.elems[i]
                    elif type(arbre.elems[j]) == Inconnue:
                        #si elles ont le même nom on les fusionne , sinon on fusionne les quantités et c'est tout
                        arbre.elems[j].quantity *= arbre.elems[i].quantity
                        arbre.elems[i].quantity = 1
                        if arbre.elems[j].name == arbre.elems[i].name:
                            arbre.elems[j].puissance += arbre.elems[i].puissance
                            del arbre.elems[i]
                        else:
                            trouve = False
                #On fusionne les nombres et les inconnues
                elif type(arbre.elems[i]) == Inconnue and type(arbre.elems[j]) == Nombre:
                    arbre.elems[i].quantity *= arbre.elems[j].value
                    del arbre.elems[j]
                    trouve = True
                j += 1
            if not trouve:
                i += 1
    return arbre

#Puissance
def compresser_puis(arbre):
    elem = arbre.elems[0]
    elem.puissance *= arbre.elems[1].value
    elem.quantity = float(elem.quantity) ** float(arbre.elems[1].value)
    return elem

#Effectuer les calculs
def operation(op,chiff1,chiff2):
    if op=='+':
        return float(chiff1)+float(chiff2)
    if op=='-':
        return float(chiff1)-float(chiff2)
    if op=='*':
        return float(chiff1)*float(chiff2)
    if op=='/':
        if float(chiff2)!=0.0:
            return float(chiff1)/float(chiff2)
    if op == "^":
        return float(chiff1)**float(chiff2)
    return "erreur"

def petits_calculs(li_toks):
    #print("on appelle petits calculs")
    li_toks = to_pnl(li_toks)
    def consume():
        operateurs = ["+","*","/","-","^"]
        if li_toks[-1] in operateurs:
            cur = li_toks.pop(-1)
            b = consume()
            a = consume()

            return operation(cur,a,b)
        cur = li_toks.pop(-1)
        return cur
    return consume()

###print("Calc1",petits_calculs( ['(',3.0,'+',4.0,'*',2.0,')',"^","(",2,")"]))
###print("Calc2", petits_calculs(["(",2,")"]))

def enlever_0(liste_to):
    #print("on appelle enlever 0")
    i=0
    if type(liste_to) == str:
        liste_to = get_tokens(liste_to)
    li_zeros = ["0.0","-0.0",0.0,-0.0]
    new_liste=[]
    taille_liste=len(liste_to)
    while i < taille_liste:

        if liste_to[i]=="(" or liste_to[i]==")":
            i+=1

        #cas avec une addition de 0
        elif liste_to[i] in li_zeros and i+1 <taille_liste and liste_to[i+1]=='+':
            i+=2

        #cas multiplication de 0
        elif liste_to[i] in li_zeros and i+1 <taille_liste and liste_to[i+1]=='*':
            i+=1
            while i < taille_liste and (liste_to[i]!='+' and liste_to[i]!='-'):
                i+=1
            i+=1

        else:
            new_liste.append(liste_to[i])
            i+=1

    liste_symb=['+','*','-']
    if new_liste[-1] in liste_symb:
        new_liste.pop()
    return new_liste

###Gestion des parenthèses
##def enlever_parentheses(chaine):
##    #print("on appelle enlever parenthese")
##    new=""
##    for lettre in chaine:
##        if lettre!="(" and lettre!=")":
##            new = new+lettre
##    return new

def trouver_l_ouvrante(token,i):
    #print("on appelle trouver l'ouvrante")
    parenthese_fermante = 1
    terme=[]
    while parenthese_fermante!=0:
        i-=1
        if token[i]==")":
            parenthese_fermante+=1
        elif token[i]=="(":
            parenthese_fermante-=1
        terme.append(token[i])
    terme.pop(-1)
    terme = terme[::-1]
    return terme,i #i position du debut du premier terme utilise pour surppimer

def trouver_la_fermante(token,i):
    #print("on appelle trouver la fermante")
    parenthese_ouvrante = 1
    terme = []
    while parenthese_ouvrante != 0:
        i += 1
        if token[i] == "(":
            parenthese_ouvrante += 1
        elif token[i] == ")":
            parenthese_ouvrante -= 1
        terme.append(token[i])
    terme.pop()
    return terme,i

#sert a trouver le terme sans parenthese par exemple 3*x1 ou 3
def terme_chiffre_droite(token,i):
    #print("on appelle terme chiffre droite")
    ok=1
    terme=[token[i]]
    while ok==1:
        i+=1
        if i < len(token) and token[i] =="*":
            terme.append(token[i])
            i+=1
            terme.append(token[i])
        else:
            ok=0
    return terme,i-1

def terme_chiffre_gauche(token,i):
    #print("on appelle terme chiffre gauche")
    ok = 1
    terme = [token[i]]
    while ok == 1:
        i -= 1
        if i > 0 and token[i] == "*":
            terme.insert(0,token[i])
            i -= 1
            terme.insert(0,token[i])
        else:
            ok = 0
    return terme,i+1


# =============================================================================
#  REPRESENTATION DES POLYNOMES
# =============================================================================

# =============================================================================
# Représentation des polynômes en tokens
# =============================================================================

## BUT DE LA PARTIE: Ecrire un polynôme (de type str) en une liste de tokens
###ordonnés en notation polonaise inversée

def get_tokens(chaine):
    #print("on appelle get tokens")
    """Transforme la chaîne de caractères "chaine" en une liste de tokens
    Exemple : "aaa+2" = ["aaa",'+',2.0] (les chiffres sont transformés en float)"""
    i = 0
    tokens = []
    operateurs = ["+","*","-","/","(",")","^"]
    #Tant qu'il nous reste un morceau à traiter
    while i < len(chaine):
        #tant qu'on est sur du vide
        while i < len(chaine) and chaine[i] == " ":
            i += 1
        #On est sur un opérateur
        if chaine[i] in operateurs:
            tokens.append(chaine[i])
            i += 1
            continue
        #tant qu'on est sur un token
        cur = chaine[i]
        i += 1
        while i < len(chaine) and chaine[i] != " " and chaine[i] not in operateurs:
            cur += chaine[i]
            i += 1
        if (cur.replace(".","")).isdigit():
            tokens.append(float(cur))
        else:
            tokens.append(cur)
    return fusion_sub(tokens)

##print("fusion sub donne",fusion_sub(get_tokens("3+-1")))
##print(get_tokens("x1*3+x2"))

def to_pnl(tokens):
    #print("on appelle to pnl")
    """Ecrit la liste de tokens en notation polonaise inversée (change l'ordre
    d'apparition des tokens dans la liste) + les parenthèses sont supprimées"""
    output = []
    stack = []
    operateurs = ["+","*","-","/","^"]
    for token in tokens:
        if type(token) == int:
            output.append(token)
        elif token in operateurs:
            while(len(stack) > 0 and stack[-1] != "(" and (precedence(stack[-1]) > precedence(token) or (precedence(stack[-1]) == precedence(token) and left_associative(token)))):
                output.append(stack.pop(-1))
            stack.append(token)
        elif token == "(":
            stack.append(token)
        elif token == ")":
            while(stack[-1] != '('):
                output.append(stack.pop(-1))
            stack.pop(-1)
        else:
            output.append(token) #c'est une inconnue
    while(len(stack) > 0):
        output.append(stack.pop(-1))
    return output

def token_to_str(token):
    m = ""
    for c in token:
        m = m + str(c)
    return m


# =============================================================================
# Représentation des polynômes en arbres
# =============================================================================

## BUT DE LA PARTIE: Ecrire un polynôme (une liste de tokens en notation
###polonaise inversée) en un arbre

#Ici on agit récursivement sur le dernier élément , si c'est un opérateur on agit
#récursivement sur ses 2 membres sinon on retourne juste sa version pythonisée (???)
def pnl_en_arbre_binaire(li_toks):
    #print("on appelle pnl en arbre binaire")
    """Convertit une liste de tokens li_toks (en notation polonaise inversée) en
    en un arbre"""
    def consume():
        operateurs = {"+":Add,"*":Mult,"/":Div,"-":Sub,"^":Puis}
        if li_toks[-1] in operateurs:
            return operateurs[li_toks.pop(-1)](consume(),consume())
        cur = li_toks.pop(-1)
        if type(cur) == float:
            return Nombre(cur)
        return Inconnue(cur)
    return consume()

#Le but de cette fonction est de fusionner les mêmes types successifs
def arbre_binaire_en_arbre(arbre):
    #print("on appelle arbre binaire en arbre")
    if type(arbre) == Nombre or type(arbre) == Inconnue: return arbre #Rien à faire
    new_elems = []

    for elem in arbre.elems:
        elem = arbre_binaire_en_arbre(elem)
        if type(elem) == type(arbre):
            for sub_elem in elem.elems:
                new_elems.append(sub_elem)
        else:
            new_elems.append(elem)
    arbre.elems = new_elems
    return arbre

def tokens_to_tree(poly):
    #print("on appelle tokens to tree")
    notation_pnl = to_pnl(poly)
    tree_tokens = pnl_en_arbre_binaire(notation_pnl)
    tree_tokens = arbre_binaire_en_arbre(tree_tokens)
    return reduce(tree_tokens)

def str_to_tree(poly):
    #print("on appelle str to tree")
    li_tokens = get_tokens(poly)
    li_sans_puis = developper_puis2(li_tokens)
    dev_poly = trouver_ou_dev_2(li_sans_puis)
    reduced_tree = tokens_to_tree(dev_poly)
    return reorganiser(reduced_tree)


# =============================================================================
# Représentation des polynômes en dictionnaires
# =============================================================================

def tree_to_dico(root):
    r, dic = tree_to_dico2(root,{})
    return dic

def tree_to_dico2(root, dico):
    #print("on appelle tree_to_dico")
    """Convertit un polynôme (écrit sous forme d'arbre) en un dictionnaire dont
    les clés sont les inconnues

    On part du principe : que les - sont stockés devant les coeffs
    que les / sont stockés dans les coeffs
    que 2*x1*x2*3*x2 = 6*x1*x2^2"""
    if type(root) == Nombre:
        add_dans_dico('1',float(root.value),dico)    #la clé = 1 car pas d'inconnue
    elif type(root) == Inconnue:
        root.puissance = float(root.puissance)
        root.quantity = float(root.quantity)
        if root.puissance != 1.0:
            add_dans_dico(root.name+"^"+str(float(root.puissance)),float(root.quantity),dico)
        else:
            add_dans_dico(root.name,float(root.quantity),dico)
    else: #c'est un opérateur [soit + soit * normalement, tout a été développé avant (les - sont dans les coefficients)]
        if type(root) == Mult:
            coeff = 1.0
            if type(root.elems[0]) == Nombre: #c'est le coeff (il n'y en a qu'un devant)
                coeff = coeff * root.elems[0].value
                del root.elems[0]
            monome = root.to_str()
            monome = enlever_0(monome)
            monome = token_to_str(monome)
            if type(monome) == Inconnue:
                monome.puissance = float(monome.puissance)
                monome.quantity = float(monome.quantity)
            if monome == "": monome = '1' #est-ce que ça doit vraiment arriver ?
            add_dans_dico(monome,coeff,dico)
        else: #si l'opérateur est '+'
            for terme in root.elems:
                if type(terme) == Inconnue:
                    terme.puissance = float(terme.puissance)
                    terme.quantity = float(terme.quantity)
                racine, dico = tree_to_dico2(terme,dico) #on écrase l'ancien dictionnaire
    return (root,dico)

def dico_to_str(dico):
    #print("on appelle dico_to_str")
    """Transforme le polynôme écrit sous la forme d'un dictionnaire en chaîne de
    caractères"""
    poly = ""
    if '1' in dico:
        poly += str(dico['1'])
        del dico['1']
        poly += " + "
    for inconnue in dico:   #on traite chaque inconnue
        dico[inconnue] = float(dico[inconnue])
        if dico[inconnue] == 0.0:
            continue
        if dico[inconnue] != 1.0:
            poly += str(dico[inconnue]) + "*"
        poly += inconnue
        poly += " + "
    return poly[:-3]

def str_to_dico(poly):
    tokens_tries = str_to_tree(poly)
    return tree_to_dico(tokens_tries)


# =============================================================================
#  OPERATIONS SUR LES POLYNOMES
# =============================================================================

# =============================================================================
# Opérations sur les coefficients
# =============================================================================

def reduce(elem):
    #print("on appelle reduce",elem.to_str())
    """Effectue des opérations simples sur les coefficients du polynôme elem
    (écrit sous forme d'arbre) comme => développer des expressions et additionner
    ou multiplier des coefficients
    IL FAUT DEVELOPPER AVANT DE REDUIRE
    """



    if type(elem) in (Nombre,Inconnue):
        return elem
    if est_expression(elem.to_str()):
        return Nombre(petits_calculs(get_tokens(elem.to_str())))
    #est_primitif = True
    for x in range(len(elem.elems)):
        elem.elems[x] = reduce(elem.elems[x])
        #if type(elem.elems[x]) != Nombre:
        #    est_primitif = False
    if type(elem) == Mult:
        elem = compresser_mult(elem)
    elif type(elem) == Add:
        elem = compresser_add(elem)
    elif type(elem) == Puis and type(elem.elems[0]) == Inconnue:
        elem = compresser_puis(elem)
    #elif est_primitif:
    #    return Nombre(petits_calculs(get_tokens(elem.to_str())))
    #print("je sors de reduce",elem.to_str())
    return elem

#Pour réorganiser un polynôme , une addition de multiplications
def reorganiser(poly):
    #print("on appelle reorganiser")
    #on assume que le type de poly == Add, on peut mettre un assert
    if type(poly) == Inconnue or type(poly) == Nombre: return poly
    i = 0
    if type(poly) == Mult:
        mult = True
        cur = poly
        poly = Add(cur)
    else:
        mult = False
    for elem in poly.elems:
        if type(elem) == Mult:
            sub_elems = []
            new_elems = []
            for x in range(len(elem.elems)):
                if elem.elems[x].quantity != 1:
                    new_elems.append(Nombre(elem.elems[x].quantity))
                    elem.elems[x].quantity = 1
                sub_elems.append((-elem.elems[x].puissance,hash(elem.elems[x].name),x))
            sub_elems.sort()
            for k in range(len(sub_elems)):
                new_elems.append(elem.elems[sub_elems[k][2]])
            elem.elems = new_elems
        i += 1
    if mult:
        poly = poly.elems[0]
    return poly


# =============================================================================
# Développer les puissances
# =============================================================================

def developper_puis2(poly):
    #print("on appelle developper puis2")
    """poly (liste de tokens) => str sans puiss"""
    poly_str = ""
    i=0
    prochain = i
    besoin_de_changer = True
    while (i <len(poly)):
        if poly[i] == '^':

            #TROUVER LA PUISSANCE
            if type(poly[i+1]) != str: #c'est un chiffre, float ou int
                a_la_puissance = float(poly[i+1])
                prochain = i+2
            elif poly[i+1] == '(': #c'est une parenthèse (
                stack_brackets = []
                j = i+1
                li_calcul = [poly[j]]
                j +=1
                while not((poly[j] == ')' and stack_brackets == [])):
                    if poly[j] == '(':
                        stack_brackets.append(poly[j])
                    elif poly[j] == ')':
                        stack_brackets.pop(-1)
                    li_calcul.append(poly[j])
                    j+=1
                li_calcul.append(poly[j]) #on ajoute la dernière parenthèse fermante )
                prochain = j + 1
                a_la_puissance = float(petits_calculs(li_calcul))
            #sinon c'est une inconnue et c'est pas normal, pas de x1^x2

            #ON TROUVE LE TERME (même si c'est un chiffre)
            if type(poly[i-1]) != str: #c'est un chiffre, float ou int
                terme = str(poly[i-1])
            elif poly[i-1] == ')': #c'est une parenthèse )
                k = i-1
                terme = ""
                while poly[k] != '(':
                    terme = str(poly[k]) + terme
                    k-=1
                terme = poly[k] +terme #on ajoute la parenthèse ouvrante (
            else:
                #sinon c'est une inconnue et on ne fait rien, x1^2 est supporté par reduce
                besoin_de_changer = False

            #ON RECOPIE n FOIS LE TERME (même si c'est un chiffre)
            if a_la_puissance > 0 and besoin_de_changer:
                p = 1
                #logiquement le terme a déjà été écrit une première fois car rencontré avant de croiser le ^
                while (p < a_la_puissance):
                    poly_str = poly_str + '*' + terme
                    p+=1
            elif not(besoin_de_changer):
                poly_str = poly_str + '^' + str(float(a_la_puissance))
            elif a_la_puissance == 0:
                for caractere in terme: #on efface la chaine qu'on a écrite une première fois
                    poly_str = poly_str[:-1]
                poly_str = poly_str + "1.0"
        else:
            poly_str += str(poly[i])
        i+=1
        i = max(prochain,i)
        besoin_de_changer = True
    return poly_str

##p = ['x1','^',2,'+','x2','^',2,'+',3,'-',2,'-',1,'+',2,'^',2,'+','(',2,'*','x3',')','^',25]
###print(developper_puis2(p))

#ne prend que deux termes
def developper_termes(arbre):
    #print("on appelle developper termes")
    if type(arbre.elems[0]) == Inconnue or type(arbre.elems[0]) == Nombre or type(arbre.elems[0]) == Mult:
        fils_1 = [arbre.elems[0]]
    else:
        fils_1 = arbre.elems[0].elems  # liste des valeurs dans le second add

    if type(arbre.elems[1])==Inconnue or type(arbre.elems[1])==Nombre or type(arbre.elems[1]) == Mult:
        fils_2=[arbre.elems[1]]
    else:
        fils_2 = arbre.elems[1].elems #liste des valeurs dans le second add
        #a forcement deux fils donc je met les deux dans des variables(qui sont forcément sous forme de liste)
    new_fils=Nombre(0)
        #new_fils va stocker la valeur du developpement ensuite je fais un double parcours pour additionner

    for pre in fils_1:
        pre_s=copy.deepcopy(pre)
        for sec in fils_2:
            sec_s=copy.deepcopy(sec)
            new_fils=Add(new_fils,Mult(pre_s,sec_s))
            pre_s = copy.copy(pre)

    return new_fils

riemann = "((x1+x2)*((x3+x4)))"
riri = get_tokens(riemann)
riri = to_pnl(riri)
riri = pnl_en_arbre_binaire(riri)
riri = arbre_binaire_en_arbre(riri)
#riri = get_tokens(riri)
riri = developper_termes(riri)
riri = get_tokens(riri.to_str())
##print(riri)
##print(riri, type(riri), riri.to_str())
new_terme = enlever_0(riri)
new_terme = token_to_str(new_terme)
#print("Riemann", new_terme)

def trouver_ou_dev_2(str_poly):
    #print("on appelle trouver ou dev 2")
    """renvoie un token"""
    token=get_tokens(str_poly)
    return trouver_ou_dev_3(token)

def trouver_ou_dev_3(token):
    #print("on appelle trouver ou dev 3")
    i=0
    while i < len(token):
        if token[i]=='*':

            if token[i+1]=="(" or token[i-1]==")":
                #Trouver les deux termes à multiplier
                if token[i-1]==")":
                    terme1,debut_terme1 = trouver_l_ouvrante(token,i-1)
                else:
                    terme1,debut_terme1 = terme_chiffre_gauche(token,i-1)
                if token[i+1]=="(":
                    terme2, fin_terme2 = trouver_la_fermante(token,i+1)
                else:
                    terme2,fin_terme2 = terme_chiffre_droite(token,i+1)

                terme1 = token[debut_terme1:i]
                terme2 = token[i+1:fin_terme2+1]

                #Tant qu'on doit développer, on agit sur les termes fils
                while cote_a_cote('*',')','(',terme1):
                    terme1 = trouver_ou_dev_3(terme1)

                while cote_a_cote('*',')','(',terme2):
                    terme2 = trouver_ou_dev_3(terme2)

                #transforme les deux termes en arbres non binaires puis les multiplie
                terme1=tokens_to_tree(terme1)
                terme2=tokens_to_tree(terme2)

                miaou_arbre = Mult(terme1,terme2)
                new_terme=developper_termes(miaou_arbre)
                new_terme = enlever_0(new_terme.to_str())
                token = token[:debut_terme1]+['(']+new_terme+[')']+token[fin_terme2+1:]

            #sinon on ne fait rien, pas de développement à effectuer
        i+=1
    return token

#paul = "-2.0*(((x1+x2)*((x3+x4))*(x5)))+3*x2"
##print("Paul", trouver_ou_dev_2(paul))


# =============================================================================
# Traiter le polynôme
# =============================================================================

def developper(poly):
    #print("on appelle developper")
    """Effectue des calculs sur la chaine de caractères poly pour faire disparaître
    les multiplications et les puissances
    Sortie : un poly prêt à être traité par Waring (str)"""
    dico = str_to_dico(poly)
    return dico_to_str(dico)

# =============================================================================
# Algorithme de Waring
# =============================================================================

def Waring(poly):
    """Applique l'algorithme de Waring à un polynôme (de type str) simplifié (sous
    forme développée)
    Renvoie le polynôme exprimé en fonction des polynômes symétriques élémentaires"""
    p, li = Waring2(poly,{}) #[] = liste des termes (du type sigma) du polynôme
    return li #à retransformer en poly

def Waring2(poly,dico_sigmas):
    #print("on appelle waring")
    p = developper(poly)
    if p == "0.0" or p =="" : return p
    li_inconnues = noms_inconnues(p)
    i = 1
    nb_inconnues = len(li_inconnues)

    #Monôme dominant
    dominant = terme_dominant(p)
    coeff_dominant = float(coeff_de(dominant, p))

    #Liste des puissances : a
    tokens_tries = str_to_tree(dominant)
    a = list_pow(tokens_tries,nb_inconnues) #liste des puissances de chaque inconnue de dominant(dans leur ordre d'apparition)

    if coeff_dominant == 1.0:
        p = p + "- (("
    else:
        p = p + "- (" + str(coeff_dominant) + '*' + '('
    affichage = "" #pour afficher sigma1^2 + ...
    puis = float(a[i-1]-a[i])
    sigmas = '(' + sigma(nb_inconnues,i) + ")^(" + str(puis) + ')'
    if puis > 0:
        affichage += "sigma" + str(i)
        if puis != 1.0:
            affichage += "^(" + str(puis) + ')'
    operateurs = ['+','-','*','/','(']
    while (i < nb_inconnues):
        i += 1
        puis = float(a[i-1]-a[i])
        if puis > 0:
            sigmas = sigmas + '*' + '('+ sigma(nb_inconnues,i) + ")" + "^(" + str(puis) + ')'
            if (affichage != "" and affichage[-1] not in operateurs):
                affichage += "*"
            affichage += "sigma" + str(i)
            if puis != 1.0:
                affichage += "^(" + str(puis) + ')'
    p = p + sigmas + ')' + ')'

    add_dans_dico(affichage,coeff_dominant,dico_sigmas) #le coeff dominant n'apparaît pas dans les sigmas
    new_p = Waring2(p, dico_sigmas)
    return (new_p, dico_sigmas)


# =============================================================================
#  Navigation
# =============================================================================

def menu_principal():
    print("Menu :")
    print("0 : Exécuter l'ensemble des tests")
    print("1 : Exécuter le programme principal")
    print("2 : Exécuter le programme principal")
    print("3 : Exécuter le programme principal")
    print("4 : Exécuter le programme principal")
    print("")

def choix_exec(msg,possibilites):
    """renvoie un choix d'exécution valide"""
    partie=input(msg)
    while partie not in possibilites:
        partie=input(msg)
    return partie


# =============================================================================
# Programme principal
# =============================================================================

menu_principal()

choix = choix_exec("Votre choix (0, 1, 2, 3 ou 4) : ",('0','1','2','3','4'))

if choix == '0': #Exécution de l'ensemble des tests
    print("Début des tests...")
    ##exec_test()
    print("Tous les tests ont été correctement effectués")
else:
    import time
    p = "(x1-x2)^2*(x1-x3)^2*(x1-x4)^2*(x2-x3)^2*(x2-x4)^2*(x3-x4)^2"
    print("p =",p)
    t1 = time.time()
    waring = Waring(p)
    t2 = time.time()
    print("Waring exécuté en",t2-t1,"secondes")
    print(waring)
    print(jolis_sigmas(dico_to_str(waring)))