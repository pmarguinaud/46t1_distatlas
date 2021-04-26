
include Makefile.inc

all: pgd

src/local/atlas/programs/atlas-pgd-demo.o: src/local/atlas/programs/atlas-pgd-demo.F90 src/local/ifsaux/module/parkind1.o src/local/atlas/io/atlas-io-gathscat.o src/local/atlas/io/atlas-io-dh.o src/local/atlas/io/atlas-fmt-null.o src/local/atlas/programs/atlas-helper.o src/local/ifsaux/module/xrd_getoptions.o
	cd src/local/atlas/programs && $(F90) -c -I../include -I../../ifsaux/module -I../io -I../io -I../io -I. -I../../ifsaux/module -o atlas-pgd-demo.o atlas-pgd-demo.F90

src/local/atlas/programs/atlas-helper.o: src/local/atlas/programs/atlas-helper.F90 src/local/ifsaux/module/parkind1.o
	cd src/local/atlas/programs && $(F90) -c -I../include -I../../ifsaux/module -o atlas-helper.o atlas-helper.F90

src/local/atlas/interpolation/interpolation4.o: src/local/atlas/interpolation/interpolation4.cc 
	cd src/local/atlas/interpolation && $(CXX) -c -I../include -o interpolation4.o interpolation4.cc

src/local/atlas/interpolation/interpolationA.o: src/local/atlas/interpolation/interpolationA.cc 
	cd src/local/atlas/interpolation && $(CXX) -c -I../include -o interpolationA.o interpolationA.cc

src/local/atlas/fortran/interpolationA_mod.o: src/local/atlas/fortran/interpolationA_mod.F90 
	cd src/local/atlas/fortran && $(F90) -c  -o interpolationA_mod.o interpolationA_mod.F90

src/local/atlas/fortran/interpolation4_mod.o: src/local/atlas/fortran/interpolation4_mod.F90 
	cd src/local/atlas/fortran && $(F90) -c  -o interpolation4_mod.o interpolation4_mod.F90

src/local/atlas/fortran/gradient_mod.o: src/local/atlas/fortran/gradient_mod.F90 
	cd src/local/atlas/fortran && $(F90) -c  -o gradient_mod.o gradient_mod.F90

src/local/atlas/gradient/gradient.o: src/local/atlas/gradient/gradient.cc 
	cd src/local/atlas/gradient && $(CXX) -c -I../include -o gradient.o gradient.cc

src/local/atlas/surfex/atlas-compute-aos.o: src/local/atlas/surfex/atlas-compute-aos.F90 src/local/ifsaux/module/parkind1.o src/local/atlas/fortran/interpolationA_mod.o src/local/atlas/fortran/interpolation4_mod.o src/local/atlas/fortran/gradient_mod.o src/local/atlas/programs/atlas-helper.o
	cd src/local/atlas/surfex && $(F90) -c -I../../ifsaux/module -I../fortran -I../fortran -I../fortran -I../programs -o atlas-compute-aos.o atlas-compute-aos.F90

src/local/atlas/surfex/atlas-compute-sso.o: src/local/atlas/surfex/atlas-compute-sso.F90 src/local/ifsaux/module/parkind1.o src/local/atlas/fortran/interpolationA_mod.o src/local/atlas/fortran/gradient_mod.o src/local/atlas/programs/atlas-helper.o
	cd src/local/atlas/surfex && $(F90) -c -I../../ifsaux/module -I../fortran -I../fortran -I../programs -o atlas-compute-sso.o atlas-compute-sso.F90

src/local/atlas/surfex/atlas-prep-impl.o: src/local/atlas/surfex/atlas-prep-impl.F90 src/local/ifsaux/module/parkind1.o src/local/atlas/fortran/interpolationA_mod.o src/local/atlas/fortran/interpolation4_mod.o src/local/atlas/io/atlas-io.o src/local/atlas/programs/atlas-helper.o
	cd src/local/atlas/surfex && $(F90) -c -I../include -I../../ifsaux/module -I../fortran -I../fortran -I../io -I../programs -o atlas-prep-impl.o atlas-prep-impl.F90

src/local/atlas/surfex/atlas-pgd-impl.o: src/local/atlas/surfex/atlas-pgd-impl.F90 src/local/ifsaux/module/parkind1.o src/local/atlas/fortran/interpolationA_mod.o src/local/atlas/fortran/interpolation4_mod.o src/local/atlas/fortran/gradient_mod.o src/local/atlas/io/atlas-io.o src/local/atlas/programs/atlas-helper.o src/local/atlas/surfex/readcovers_mod.o
	cd src/local/atlas/surfex && $(F90) -c -I../include -I../../ifsaux/module -I../fortran -I../fortran -I../fortran -I../io -I../programs -I. -o atlas-pgd-impl.o atlas-pgd-impl.F90

src/local/atlas/surfex/readcovers_mod.o: src/local/atlas/surfex/readcovers_mod.F90 src/local/ifsaux/module/parkind1.o
	cd src/local/atlas/surfex && $(F90) -c -I../../ifsaux/module -o readcovers_mod.o readcovers_mod.F90

src/local/atlas/surfex/atlas-compute-covers.o: src/local/atlas/surfex/atlas-compute-covers.F90 src/local/ifsaux/module/parkind1.o src/local/atlas/fortran/interpolationA_mod.o src/local/atlas/programs/atlas-helper.o
	cd src/local/atlas/surfex && $(F90) -c -I../include -I../../ifsaux/module -I../fortran -I../programs -o atlas-compute-covers.o atlas-compute-covers.F90

src/local/atlas/io/atlas-fmt-null.o: src/local/atlas/io/atlas-fmt-null.F90 src/local/atlas/io/atlas-fmt.o src/local/ifsaux/module/parkind1.o
	cd src/local/atlas/io && $(F90) -c -I../include -I. -I../../ifsaux/module -o atlas-fmt-null.o atlas-fmt-null.F90

src/local/atlas/io/atlas-io-gathscat.o: src/local/atlas/io/atlas-io-gathscat.F90 src/local/atlas/io/atlas-io.o src/local/atlas/io/atlas-fmt.o src/local/ifsaux/module/parkind1.o
	cd src/local/atlas/io && $(F90) -c -I../include -I. -I. -I../../ifsaux/module -o atlas-io-gathscat.o atlas-io-gathscat.F90

src/local/atlas/io/atlas-io-dh.o: src/local/atlas/io/atlas-io-dh.F90 src/local/atlas/io/atlas-io.o src/local/ifsaux/module/parkind1.o
	cd src/local/atlas/io && $(F90) -c -I../include -I. -I../../ifsaux/module -o atlas-io-dh.o atlas-io-dh.F90

src/local/atlas/io/atlas-fmt.o: src/local/atlas/io/atlas-fmt.F90 src/local/ifsaux/module/parkind1.o
	cd src/local/atlas/io && $(F90) -c -I../../ifsaux/module -o atlas-fmt.o atlas-fmt.F90

src/local/atlas/io/atlas-io.o: src/local/atlas/io/atlas-io.F90 src/local/ifsaux/module/parkind1.o
	cd src/local/atlas/io && $(F90) -c -I../../ifsaux/module -o atlas-io.o atlas-io.F90

src/local/ifsaux/module/parkind1.o: src/local/ifsaux/module/parkind1.F90 
	cd src/local/ifsaux/module && $(F90) -c  -o parkind1.o parkind1.F90

src/local/ifsaux/module/xrd_getoptions.o: src/local/ifsaux/module/xrd_getoptions.F90 src/local/ifsaux/module/parkind1.o src/local/ifsaux/module/xrd_unix_env.o
	cd src/local/ifsaux/module && $(F90) -c -I. -I. -o xrd_getoptions.o xrd_getoptions.F90

src/local/ifsaux/module/xrd_unix_env.o: src/local/ifsaux/module/xrd_unix_env.F90 src/local/ifsaux/module/parkind1.o
	cd src/local/ifsaux/module && $(F90) -c -I. -o xrd_unix_env.o xrd_unix_env.F90

src/local/ifsaux/fi_libc/fi_libc.o: src/local/ifsaux/fi_libc/fi_libc.c 
	cd src/local/ifsaux/fi_libc && $(CC) -c -I. -o fi_libc.o fi_libc.c

src/local/ifsaux/support/iswap8.o: src/local/ifsaux/support/iswap8.c 
	cd src/local/ifsaux/support && $(CC) -c  -o iswap8.o iswap8.c


pgd: src/local/atlas/programs/atlas-pgd-demo.o src/local/atlas/programs/atlas-helper.o src/local/atlas/interpolation/interpolation4.o src/local/atlas/interpolation/interpolationA.o src/local/atlas/fortran/interpolationA_mod.o src/local/atlas/fortran/interpolation4_mod.o src/local/atlas/fortran/gradient_mod.o src/local/atlas/gradient/gradient.o src/local/atlas/surfex/atlas-compute-aos.o src/local/atlas/surfex/atlas-compute-sso.o src/local/atlas/surfex/atlas-prep-impl.o src/local/atlas/surfex/atlas-pgd-impl.o src/local/atlas/surfex/readcovers_mod.o src/local/atlas/surfex/atlas-compute-covers.o src/local/atlas/io/atlas-fmt-null.o src/local/atlas/io/atlas-io-gathscat.o src/local/atlas/io/atlas-io-dh.o src/local/atlas/io/atlas-fmt.o src/local/atlas/io/atlas-io.o src/local/ifsaux/module/parkind1.o src/local/ifsaux/module/xrd_getoptions.o src/local/ifsaux/module/xrd_unix_env.o src/local/ifsaux/fi_libc/fi_libc.o src/local/ifsaux/support/iswap8.o
	$(LD) -o pgd src/local/atlas/programs/atlas-pgd-demo.o src/local/atlas/programs/atlas-helper.o src/local/atlas/interpolation/interpolation4.o src/local/atlas/interpolation/interpolationA.o src/local/atlas/fortran/interpolationA_mod.o src/local/atlas/fortran/interpolation4_mod.o src/local/atlas/fortran/gradient_mod.o src/local/atlas/gradient/gradient.o src/local/atlas/surfex/atlas-compute-aos.o src/local/atlas/surfex/atlas-compute-sso.o src/local/atlas/surfex/atlas-prep-impl.o src/local/atlas/surfex/atlas-pgd-impl.o src/local/atlas/surfex/readcovers_mod.o src/local/atlas/surfex/atlas-compute-covers.o src/local/atlas/io/atlas-fmt-null.o src/local/atlas/io/atlas-io-gathscat.o src/local/atlas/io/atlas-io-dh.o src/local/atlas/io/atlas-fmt.o src/local/atlas/io/atlas-io.o src/local/ifsaux/module/parkind1.o src/local/ifsaux/module/xrd_getoptions.o src/local/ifsaux/module/xrd_unix_env.o src/local/ifsaux/fi_libc/fi_libc.o src/local/ifsaux/support/iswap8.o $(LIBS)

clean: 
	\rm -f src/local/atlas/programs/atlas-pgd-demo.o src/local/atlas/programs/atlas-helper.o src/local/atlas/interpolation/interpolation4.o src/local/atlas/interpolation/interpolationA.o src/local/atlas/fortran/interpolationA_mod.o src/local/atlas/fortran/interpolation4_mod.o src/local/atlas/fortran/gradient_mod.o src/local/atlas/gradient/gradient.o src/local/atlas/surfex/atlas-compute-aos.o src/local/atlas/surfex/atlas-compute-sso.o src/local/atlas/surfex/atlas-prep-impl.o src/local/atlas/surfex/atlas-pgd-impl.o src/local/atlas/surfex/readcovers_mod.o src/local/atlas/surfex/atlas-compute-covers.o src/local/atlas/io/atlas-fmt-null.o src/local/atlas/io/atlas-io-gathscat.o src/local/atlas/io/atlas-io-dh.o src/local/atlas/io/atlas-fmt.o src/local/atlas/io/atlas-io.o src/local/ifsaux/module/parkind1.o src/local/ifsaux/module/xrd_getoptions.o src/local/ifsaux/module/xrd_unix_env.o src/local/ifsaux/fi_libc/fi_libc.o src/local/ifsaux/support/iswap8.o pgd


