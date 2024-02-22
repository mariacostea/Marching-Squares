# Marching-Squares
Pentru a putea paraleliza funcțiile am creat o structura cu parametrii funcțiilor
pe care urma să le paralelizez. În main am inițializat câmpurile structurii, apoi 
am creat threadurile. Am paralelizat cele trei funcții: rescale, sample_grid și 
march în ordinea în care erau apelate la final în main, deoarece fiecare depindea 
de rezultatul funcției anterioare. Am calculat pentru fiecare funcție paralelizată
start-ul și end-ul porțiunii din imagine pe care a trebuit să o modifice fiecare 
thread dat ca argument. Bariera așteaptă ca toate threadurile să termine de 
redimensionat partea lor asignată din imagine pentru a putea continua cu modificările 
aduse de următoarele funcții. Apoi se paralelizează funcția care asignează valori 
binare matricii și funcția care face desenele în funcție de acele valori binare din 
colțuri, din nou fiecare thread este responsabil să aducă aceste modificări părții 
din imaginea finală, care este delimitată de start-ul și end-ul calculat pentru 
funcția respectivă. La final, valorile care erau returnate de funcțiile paralelizate
sunt salvate în structură și am modificat ultimele linii care apelau aceste funcții.
