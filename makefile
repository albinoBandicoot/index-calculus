indcalc: *.c *.h
	gcc -O3 -o indcalc tdiv.c relation.c matsolve.c indcalc.c core.c -I ./ -lgmp

safegen: safegen.c
	gcc -o safegen safegen.c -lgmp

clean:
	rm -f indcalc
