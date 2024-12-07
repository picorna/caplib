.SUFFIXES:
.SUFFIXES: .c .o
.IGNORE:

.c.o:
	$(CC) -I$(INC) $(AUX_FLAGS) $<

CC = cc -c -DDEBUG=0 -O0
LINK = cc -lm -O0
AUX_FLAGS =
INC  = ./include
CFLAG = -I$(INC)
LFLAG =
EXE = ../bin

# Selection:
all:
	@echo 'make capaxis'
	@echo 'make capgen'
	@echo 'make cappack'
	@echo 'make caprot'
	@echo 'make clean     - .o files are removed'

capgen:
	make load_capgen PNAME="capgen" AUX_FLAGS=""

capaxis:
	make load_capaxis PNAME="capaxis" AUX_FLAGS=""

cappack:
	make load_cappack PNAME="cappack" AUX_FLAGS=""

caprot:
	make load_caprot PNAME="caprot" AUX_FLAGS=""

# Clean all the objects
clean:
	rm -f $(OBJS_CAPGEN)
	rm -f $(OBJS_CAPROT)
	rm -f $(OBJS_CAPAXIS)
	rm -f $(OBJS_CAPPACK)

EXE  = ../bin

OBJS_CAPGEN = capgen_main.o capgen_call_getopt.o  capgen_gen_entire.o \
	analyse_biomt.o calc_cell.o \
	cell_axis_fit.o cell_axis_fit_any.o \
	cell_number.o datashift.o eqsolve.o error_print.o error_store.o \
	for_a_pdb.o front.o \
	make_orthogonal_matrix.o matrix_determinant.o matrix_inverse.o \
	matrix_product.o print_a_pdb.o print_center.o print_rot.o \
	read_a_pdb.o read_biomt_file.o read_wat216.o rotate_atoms.o \
	set_2plv_rotmat.o set_bplane.o set_planes.o set_ptable.o \
	sort_double.o summary_print.o summary_store.o trans_2plvto4rhv.o capgen_usage.o

OBJS_CAPAXIS = capaxis_main.o capaxis_call_getopt.o add_wat216_to_2PLV.o analyse_biomt.o calc_cell.o \
	cell_axis_fit.o cell_axis_fit_any.o \
	cell_number.o datashift.o eqsolve.o error_print.o error_store.o \
	for_a_pdb.o front.o gen_entire.o \
	make_orthogonal_matrix.o matrix_determinant.o matrix_inverse.o \
	matrix_product.o print_a_pdb.o print_center.o print_rot.o \
	read_a_pdb.o read_biomt_file.o read_wat216.o rotate_atoms.o \
	set_2plv_rotmat.o set_bplane.o set_planes.o set_ptable.o \
	sort_double.o summary_print.o summary_store.o trans_2plvto4rhv.o usage.o

OBJS_CAPPACK = cappack_main.o call_cappack.o cappack_call_getopt.o  capgen_gen_entire.o \
	analyse_biomt.o calc_cell.o \
	cell_axis_fit.o cell_axis_fit_any.o \
	cell_number.o datashift.o eqsolve.o error_print.o error_store.o \
	for_a_pdb.o front.o \
	make_orthogonal_matrix.o matrix_determinant.o matrix_inverse.o \
	matrix_product.o print_a_pdb.o print_center.o print_rot.o \
	read_a_pdb.o read_biomt_file.o read_wat216.o rotate_atoms.o \
	set_2plv_rotmat.o set_bplane.o set_planes.o set_ptable.o \
	sort_double.o summary_print.o summary_store.o trans_2plvto4rhv.o cappack_usage.o

OBJS_CAPROT = caprot_main.o caprot_call_getopt.o  gen_entire.o \
	analyse_biomt.o calc_cell.o \
	cell_axis_fit.o cell_axis_fit_any.o \
	cell_number.o datashift.o eqsolve.o error_print.o error_store.o \
	for_a_pdb.o front.o \
	make_orthogonal_matrix.o matrix_determinant.o matrix_inverse.o \
	matrix_product.o print_a_pdbs.o print_center.o print_rot.o \
	read_a_pdb.o read_biomt_file.o read_wat216.o rotate_atoms.o \
	set_2plv_rotmat.o set_bplane.o set_planes.o set_ptable.o \
	sort_double.o summary_print.o summary_store.o trans_2plvto4rhv.o caprot_usage.o

load_capaxis: $(OBJS_CAPAXIS)
	$(LINK) $(LFLAG) $(OBJS_CAPAXIS)
	mv a.out $(EXE)/$(PNAME)

load_capgen: $(OBJS_CAPGEN)
	$(LINK) $(LFLAG) $(OBJS_CAPGEN)
	mv a.out $(EXE)/$(PNAME)

load_cappack: $(OBJS_CAPPACK)
	$(LINK) $(LFLAG) $(OBJS_CAPPACK)
	mv a.out $(EXE)/$(PNAME)

load_caprot: $(OBJS_CAPROT)
	$(LINK) $(LFLAG) $(OBJS_CAPROT)
	mv a.out $(EXE)/$(PNAME)

$(OBJS_CAPGEN): Makefile
$(OBJS_CAPAXIS): Makefile
$(OBJS_CAPROT): Makefile
$(OBJS_CAPPACK): Makefile
