require "HALOBase.pm";
@myFileList=("Particles/typesDef_PIC.f90",
             "Particles/PIC/Boundary_PIC.f90",
             "Particles/PIC/ParticleInCell.f90",
             "Particles/PIC/pic_init.f90",
             "Particles/PIC/pic_tools.f90",
             "Particles/PIC/pic_treatment.f90",
             "Particles/PIC/pic_output.f90",
             "Particles/PIC/SetParticlePosition.F90");

push(@flist,@myFileList);
