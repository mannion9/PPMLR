Compiler = 'gfortran'  # Fortran compiler
COMP     = 'PC'       # Computer type

# Definitions for different systems
if COMP == 'Mac':
    Open = 'open'       # Command to open files in shell
    Mod  = '.md'       # file end for modules
    Exe  = 'a.out'      # file end for executable
    Run  = './'         # Command to run executable
if COMP == 'PC':
    Open = 'cygstart'  # Command to open files in shell
    Mod  = '.mod'       # file end for modules
    Exe  = 'a.exe'     # file end for executable
    Run  = './'         # Command to run executable

# Make   
from subprocess import call
print('  ')
print('Removing all old modules....')
a=call(['rm','*'+Mod])                   # Remove all modules for new build
print('  ')
print('Removing all old modules in RP-Exact....')
b=call(['rm','RP-Exact/*'+Mod])          # Remove all modules for new build
print('  ')
print('Removing all executables....')
c=call(['rm',Exe])                       # Remove old executible
print('  ')
print('Compiling Main program....')
d=call([Compiler,'Main.f90'])            # Run Main program
print('  ')
print('Executing....')
e=call([Run+Exe])                        # Run executable
print('  ')
print('Giving RP-Exact access to CommonData....')
f=call(['cp','*'+Mod,'RP-Exact/'])       # Allo RP-Exact to get paramters from Main.f90
print('  ')
print('Compiling exact solution....')
g=call([Compiler,'RP-Exact/Exact.f90'])  # Run Exact solution
print('  ')
print('Executing....')
h=call([Run+Exe])                        # Run executable
print('  ')
print('Plotting....')
i=call(['python','plot.py'])             # Plot data
print('  ')
print('Saving video....')
j=call([Open,'Output/out.mp4'])          # Open video

# Check For Success
success = 0
result = [c,d,e,f,g,h,i,j] # a and b dont matter just a rm statment
for i in range(len(result)):
    if result[i]!=0:
        success=1
        print('There was an issue on evaluation %i ' % (2+i)) # offshit because we ignore a , b

print('')
print('')
print('')
if success == 0:
    print('Compilation was succesful!')
else: 
    print('Compilation was not succesful.')