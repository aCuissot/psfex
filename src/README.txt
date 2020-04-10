ajouter dans le makefile généré:
.c.o:
	ajouter l'option -fPIC
(faire de même dans les repertoires fits levmar et wcs


pythonLib: setup.py libpsfexlib.a
	python3 setup.py build_ext --inplace

libpsfexlib.a:
	ar rcs libpsfexlib.a *.o fits/*.o levmar/*o wcs/*.o
