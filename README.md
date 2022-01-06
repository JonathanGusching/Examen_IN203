# Examen_machine_2022
Examen machine du 6 Janvier 2022 : compression d'image parallèle par série de Fourier

## Mesure initiale de temps 
Pour 10% de pixels et l'image tiny_lena_gray, on trouve: 
```

Transformée de Fourier, durée : 0.500444
Nombre de coefficients de fourier restant : 409
Choix des pixels, durée : 0.00039734
Reconstruction, durée : 0.070869
```

Pour 90% sur la même image :
```
Caracteristique de l'image : 64x64 => 4096 pixels.
Transformée de Fourier, durée : 0.497717
Nombre de coefficients de fourier restant : 3686
Choix des pixels, durée : 0.000987915
Reconstruction, durée : 0.562418
```

La reconstruction de l'image et le calcul de la DFT prennent ainsi beaucoup de temps, alors que la sélection des pixels est relativement négligeable.

## Parallélisation OMP
Ce qui prend le plus de temps, c'est le calcul de DFT et la reconstruction. Le choix des pixels est négligeable, on ne s'en occupe donc pas (peu de gains possibles).
Pour small_tiny_lena_gray avant parallélisation sur DFT, à 10%:
```
Transformée de Fourier, durée : 138.903
Nombre de coefficients de fourier restant : 6553
Choix des pixels, durée : 0.00695942
Reconstruction, durée : 17.039
```

Pour small_tiny_lena_gray après parallélisation sur DFT (fractionnement de la boucle qui calcule des fonctions trigonométriques, exponentielles, ...) (espace de l'image) :
```
Caracteristique de l'image : 256x256 => 65536 pixels.
Transformée de Fourier, durée : 39.9179
Nombre de coefficients de fourier restant : 6553
Choix des pixels, durée : 0.0071273
Reconstruction, durée : 16.7611
```

Idem, mais après une parallélisation sur la reconstruction (espace des fréquences):
```
Caracteristique de l'image : 256x256 => 65536 pixels.
Transformée de Fourier, durée : 39.1242
Nombre de coefficients de fourier restant : 6553
Choix des pixels, durée : 0.0069832
Reconstruction, durée : 4.9506
```

On obtient une accélération de 3.522 avec 4 threads, ce qui est très convenable.

## Parallélisation MPI 1 (espace):
Le temps de calcul est nettement réduit par rapport à la version non parallélisée, on obtient quelque chose de similaire à OpenMP, pour 4 processus:
```
Transformée de Fourier, durée : 38.3916
Transformée de Fourier, durée : 39.1456
Transformée de Fourier, durée : 39.3093
Transformée de Fourier, durée : 39.7825

Nombre de coefficients de fourier restant : 6553
Choix des pixels, durée : 0.00757006
Reconstruction, durée : 19.9268 (non parallélisée)
```

## Parallélisation MPI 2 (fréquences):
Il y a une grosse amélioration au niveau du choix des coefficients (car on a parallélisé cette partie, dans les fréquences, puis le processus 0 reconstruit l'image):
```
Caracteristique de l'image : 256x256 => 65536 pixels.
Caracteristique de l'image : 256x256 => 65536 pixels.
Caracteristique de l'image : 256x256 => 65536 pixels.
Caracteristique de l'image : 256x256 => 65536 pixels.
Transformée de Fourier, durée : 153.474
Nombre de coefficients de fourier restant : 6553
Choix des pixels, durée : 0.00312872
Transformée de Fourier, durée : 153.798
Nombre de coefficients de fourier restant : 6553
Choix des pixels, durée : 0.00288279
Transformée de Fourier, durée : 154.647
Nombre de coefficients de fourier restant : 6553
Choix des pixels, durée : 0.00269156
Transformée de Fourier, durée : 156.25
Nombre de coefficients de fourier restant : 6553
Choix des pixels, durée : 0.00285627
Reconstruction, durée : 20.3467

```

Mais cette amélioration des performances n'est pas significative en ce sens où elle ne compose qu'une petite partie du programme en terme en temps de traitement.
Avec un peu plus de temps, nous pourrions aussi éventuellement paralléliser inversePartialDiscretTransformFourier de la même manière qu'à la question précédente en espace, pour la reconstruction de l'image car celle-ci prend beaucoup de temps, ou même coupler MPI avec OpenMP pour bénéficier des deux parallélismes.
