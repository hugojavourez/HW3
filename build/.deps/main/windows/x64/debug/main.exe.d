{
    files = {
        [[build\.objs\main\windows\x64\debug\bcdata.cpp.obj]],
        [[build\.objs\main\windows\x64\debug\grid.cpp.obj]],
        [[build\.objs\main\windows\x64\debug\main.cpp.obj]],
        [[build\.objs\main\windows\x64\debug\math.cpp.obj]],
        [[build\.objs\main\windows\x64\debug\solver.cpp.obj]]
    },
    values = {
        [[C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.31.31103\bin\HostX64\x64\link.exe]],
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