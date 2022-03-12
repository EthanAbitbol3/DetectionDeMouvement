Ce projet de vision par ordinateur est divisé en 3 parties.
La première partie algorithmique demande l’implémentation d’une chaîne de traitement de type système embarqué que nous allons voir avec sigma delta.
La seconde consiste a réaliser un ensemble d’optimisations sur des morphologies
mathematiques à la fois algorithmiques, logicielles et architecturales afin d’avoir 
un rendu en temps réel et aussi de comprendre leur gain de performance séparé et
combiné.

Enfin, le dernier objectif consiste a d´evelopper une methodologie d’implantation ´
et validation pour ameliorer ses capacités de débogage.
Notre travail ici consistera tout d’abord a implémenter sigma delta suivi d’un
algorithme de morphologie “naıf“ pour obtenir une version de reférence, qu’on
utilisera dans un second temps afin de tester et d’optimiser cette chaıne de traitement
pour qu’elle soit la plus rapide possible en appliquant des transformations de bas
niveaux tel que les deroulements de boucles et des transformations de haut niveaux
: pipeline d’operateurs, fusions d’opérateurs. Pour finir, nous implémenterons du
subword parallelism (SWP) qui semble être l’optimisation la plus prometteuse.
Ainsi, nous testerons et nous comparerons toutes les fonctions pour en d´eduire la
meilleur.


Makefile: pour compiler et générer un exécutable de test
Makefile_lib: pour compiler et générer une bibliothèque qui servira pour l'évaluation

Pensez a systématiquement faire un clean

Pour mettre au point et pour compiler plus rapidement: -O0 
Pour générer un exécutable et une bibliothèque rapide: -O3

Les noms des fonctions sont ceux indiqués dans .c et .h
Ces noms sont ceux appelés par la fonction de test et la fonction de benchmark.
Ne pas changer ces noms car l'évaluation se fera pas une autre fonction morpho_test.c que celle fourni 
