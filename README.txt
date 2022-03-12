Makefile: pour compiler et générer un exécutable de test
Makefile_lib: pour compiler et générer une bibliothèque qui servira pour l'évaluation

Pensez a systématiquement faire un clean

Pour mettre au point et pour compiler plus rapidement: -O0 
Pour générer un exécutable et une bibliothèque rapide: -O3

Les noms des fonctions sont ceux indiqués dans .c et .h
Ces noms sont ceux appelés par la fonction de test et la fonction de benchmark.
Ne pas changer ces noms car l'évaluation se fera pas une autre fonction morpho_test.c que celle fourni 