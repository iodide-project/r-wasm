For this document, install ATLAS CBLAS/LAPACK. You could compile your own
CBLAS/CLAPACK.

Download F2Clib from https://www.netlib.org/f2c/. We need f2clib(the library),
make a folder and unzip this file inside the folder (else you'll clutter your workspace).

R was downloaded from the R website

Directory structure

```
Repo
|
|--> README.md
|--> R-3.6.1 unzipped
|--> libf2c/
|--> a folder called pref
```


Read README in f2clib. For hala[1], i had to run 

```
        mv f2c.h f2c.h0
        sed 's/long int /int /' f2c.h0 >f2c.h
```

(see the README file in that folder). This is important else all R calls to
fortran result in segfaults of random errors. ALSO compile with `CFLAGS=-fPIC` (by editing makefile.u and saving as Makefile)


```
libf2c/f2c.h
libf2c/libf2c.a
```

Now, go to the R folder!

In the commit (and thereafter) i have files edited to reflect changes

Essentially:

```
cp /home/sguha/mz/r-wasm/R1/R-3.6.1/share/make/vars.mk share/make/
cp /home/sguha/mz/r-wasm/R1/R-3.6.1/src/appl/dsvdc.f dtrsl.f src/appl/
cp /home/sguha/mz/r-wasm/R1/R-3.6.1/src/appl/{ dtrsl.c dpofa.c dqrdc2.c dqrutl.c dpbfa.c dqrsl.c dqrdc.c dpbsl.c dsvdc.c dtrco.c dqrls.c src/appl/
cp /home/sguha/mz/r-wasm/R1/R-3.6.1/src/appl/Makefile.in 

/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/base/baseloader.R
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/base/R/kappa.R
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/base/R/qr.R

/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/stats/R/ARMAtheory.R
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/stats/R/ar.R
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/stats/src/init.c
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/stats/src/lm.c
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/stats/src/loessc.c
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/stats/src/Makefile.in
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/library/stats/src/bsplvd.c lminfl.c kmns.c bvalue.c bvalus.c loessf.c qsbart.c sgram.c sinerp.c
sslvrg.c stxwx.c eureka.c stl.c

  
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/main/main.c
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/main/Makefile.in
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/main/registration.c
/home/sguha/mz/r-wasm/R1/R-3.6.1/src/main/Rmain.c
```

You wont need to copy anything since commit [4afa62b81308c72ed4aa72c23af88113401138a2](https://github.com/iodide-project/r-wasm/commit/4afa62b81308c72ed4aa72c23af88113401138a2) 
onward have the files modified.



Next, inside R-3.6.1, (i made a `pref` folder above this). Change the paths to reflect your installation (you need ATLAS and CBLAS )

```
CPPFLAGS="-I$HOME/mz/r-wasm/R2/libf2c" CFLAGS="-I$HOME/mz/r-wasm/R2/libf2c" MAIN_LDFLAGS="-L$HOME/mz/r-wasm/R2/libf2c/" SHLIB_LDFLAGS="-L$HOME/mz/r-wasm/R2/libf2c/" LDFLAGS="-L$HOME/mz/r-wasm/R2/libf2c/" \ 
./configure --prefix=$HOME/r-wasm/R2/pref/  --with-blas="-L/usr/lib64/atlas/ -ltatlas"     --with-lapack  --with-x=no --enable-java=no --with-readline=no --with-recommended-packages=no  --enable-BLAS-shlib=no --enable-R-shlib=yes --with-tcltk=no
make V=1 -j15 | tee -a what_run.txt
make install help

```

Confirm that `what_run.txt` does not contain any calls to `gfortran` (except for
-lgfortran which i *could not remove*).

Now you can run some tests! See https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Testing-a-Unix_002dalike-Installation
Unfortunately for some reason I've not been able to install help files fro the stats library. the tests depend on those help files being present.


1. output of /proc/cpuinfo

```
    processor       : 15
    vendor_id       : GenuineIntel
    cpu family      : 6
    model           : 85
    model name      : Intel(R) Xeon(R) Platinum 8175M CPU @ 2.50GHz
    stepping        : 4
    microcode       : 0x2000043
    cpu MHz         : 2500.000
    cache size      : 33792 KB
    physical id     : 0
    siblings        : 16
    core id         : 7
    cpu cores       : 8
    apicid          : 15
    initial apicid  : 15
    fpu             : yes
    fpu_exception   : yes
    cpuid level     : 13
    wp              : yes
    flags           : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 ss ht syscall nx pdpe1gb rdtscp lm constant_tsc rep_good nopl xtopology nonstop_tsc aperfmperf eagerfpu pni pclmulqdq ssse3 fma cx16 pcid sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand hypervisor lahf_lm abm 3dnowprefetch fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm mpx avx512f avx512dq rdseed adx smap clflushopt clwb avx512cd avx512bw avx512vl xsaveopt xsavec xgetbv1 ida arat pku ospke
    bogomips        : 5000.00
    clflush size    : 64
    cache_alignment : 64
    address sizes   : 46 bits physical, 48 bits virtual
    power management:
```
