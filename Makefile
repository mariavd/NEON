CC      = gcc
CFLAGS  = -W -Wall -Wextra

PROFILING = -pg
# note: remove optimization in order to run the profiler
# afterwards run gprof as:
#  gprof ./schrodinger > gprof.out

OPTIMIZE = -O2 -finline-functions-called-once -ftree-vectorizer-verbose=7
PEDANTIC = -pedantic -ansi 
TRAD = -Wtraditional
SHADOW = -Wshadow
ERRORS = -Werror -Wpointer-arith -Wconversion -Wcast-qual -Wwrite-strings -Wstrict-prototypes -fshort-enums -Wcast-align -fno-common -Wmissing-prototypes -Wnested-externs
IGNORE = -Wno-sign-compare -Wno-unused-result
#PUTYPE = -march=skylake
CPUTYPE =  -mavx2 -mfma -march=native 
#CPUTYPE = -march=pentium3

# -Wfatal-errors
#      This option causes the compiler to abort compilation on the first error occurred rather than trying to keep going and printing further error messages. 
# -Wunused-parameter
#  Warn whenever a function parameter is unused aside from its declaration. 
#  -Wfloat-equal
#  Warn if floating point values are used in equality comparisons. 
#  -Wundef
#  Warn if an undefined identifier is evaluated in an `#if' directive. 
# -Wunsafe-loop-optimizations
#  Warn if the loop cannot be optimized because the compiler could not assume anything on the bounds of the loop indices. With -funsafe-loop-optimizations warn if the compiler made such assumptions. 
#  -Wbad-function-cast (C only)
#      Warn whenever a function call is cast to a non-matching type. For example, warn if int malloc() is cast to anything *. # -Wstrict-prototypes (C only)
#    Warn if a function is declared or defined without specifying the argument types. (An old-style function definition is permitted without a warning if preceded by a declaration which specifies the argument types.) 
#    -Wmissing-prototypes (C only)
#        Warn if a global function is defined without a previous prototype declaration. This warning is issued even if the definition itself provides a prototype. The aim is to detect global functions that fail to be declared in header files.
#        -Wmissing-declarations (C only)
#            Warn if a global function is defined without a previous declaration. Do so even if the definition itself provides a prototype. Use this option to detect global functions that are not declared in header files. 
# -Wmissing-noreturn
#      Warn about functions which might be candidates for attribute noreturn. Note these are only possible candidates, not absolute ones. Care should be taken to manually verify functions actually do not ever return before adding the noreturn attribute, otherwise subtle code generation bugs could be introduced. You will not get a warning for main in hosted C environments.
#      -Wmissing-format-attribute
#          Warn about function pointers which might be candidates for format attributes. Note these are only possible candidates, not absolute ones. GCC will guess that function pointers with format attributes that are used in assignment, initialization, parameter passing or return statements should have a corresponding format attribute in the resulting type. I.e. the left-hand side of the assignment or initialization, the type of the parameter variable, or the return type of the containing function respectively should also have a format attribute to avoid the warning. 
#  -Wpacked
#      Warn if a structure is given the packed attribute, but the packed attribute has no effect on the layout or size of the structure. Such structures may be mis-aligned for little benefit. 
#  -Wpadded
#      Warn if padding is included in a structure, either to align an element of the structure or to align the whole structure. Sometimes when this happens it is possible to rearrange the fields of the structure to reduce the padding and so make the structure smaller.
#      -Wredundant-decls
#          Warn if anything is declared more than once in the same scope, even in cases where multiple declaration is valid and changes nothing.
#          -Wnested-externs (C only)
#              Warn if an extern declaration is encountered within a function.
#              -Wunreachable-code
#                  Warn if the compiler detects that code will never be executed. 
# -Winline
#      Warn if a function can not be inlined and it was declared as inline. Even with this option, the compiler will not warn about failures to inline functions declared in system headers.
#
#          The compiler uses a variety of heuristics to determine whether or not to inline a function. For example, the compiler takes into account the size of the function being inlined and the amount of inlining that has already been done in the current function. Therefore, seemingly insignificant changes in the source program can cause the warnings produced by -Winline to appear or disappear. 
#

#ERRORS += $(PEDANTIC)
#ERRORS += $(TRAD)
#ERRORS += $(SHADOW)

#CFLAGS += $(ERRORS)

#CFLAGS += $(PROFILING)
CFLAGS += $(OPTIMIZE)
CFLAGS += $(CPUTYPE)
CFLAGS += $(IGNORE)
CFLAGS += $(VERSION)

TARGET = schrodinger

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) $< -o $@ -lm -lgsl -lgslcblas

clean:
	rm -f $(TARGET)

