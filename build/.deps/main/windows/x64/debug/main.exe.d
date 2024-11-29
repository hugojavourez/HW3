{
    files = {
        [[build\.objs\main\windows\x64\debug\grid.cpp.obj]],
        [[build\.objs\main\windows\x64\debug\main.cpp.obj]],
        [[build\.objs\main\windows\x64\debug\manager.cpp.obj]],
        [[build\.objs\main\windows\x64\debug\math.cpp.obj]],
        [[build\.objs\main\windows\x64\debug\solution.cpp.obj]],
        [[build\.objs\main\windows\x64\debug\solver.cpp.obj]]
    },
    values = {
        [[C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Tools\MSVC\14.42.34433\bin\HostX64\x64\link.exe]],
        {
            "-nologo",
            "-dynamicbase",
            "-nxcompat",
            "-machine:x64",
            "-debug",
            [[-pdb:build\windows\x64\debug\main.pdb]]
        }
    }
}