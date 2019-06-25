FC = gfortran

FF = -g -O3 -fbounds-check -fbacktrace

# Path to the netCDF libraries
NCDF_ROOT = /usr/local/Cellar/netcdf/4.6.2
LIB = -L$(NCDF_ROOT)/lib -lnetcdf -lnetcdff

# Name of the file to compile
IZEN = carbon_flux
#IZEN = carbon_accu
#IZEN = carbon_psi
#IZEN = carbon_vol

NAME = $(IZEN).f90

All: $(IZEN).x

$(IZEN).x: $(NAME)
	$(FC) $(FF)  -o $(IZEN).x $(NAME) -I $(NCDF_ROOT)/include $(LIB)

hej:
	@echo hej!


clean:
	rm -f *.x #

