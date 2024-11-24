set_project("__name__")
set_xmakever("2.8.6")

add_rules("mode.debug", "mode.release", "mode.releasedbg")

set_policy("build.sanitizer.address", false)
set_policy("build.warning", true)
set_policy("run.autobuild", true)
set_policy("build.optimization.lto", false)

set_warnings("all") 
set_languages("c++17", "c99") 
set_runtimes("MD") -- msvc runtime library (MD/MT/MDd/MTd)

add_requires("eigen")

target("main")
    set_kind("binary")
    add_packages("eigen")
    add_files("*.cpp")
   
-- xmake project -k compile_commands  -- pour ajouter une librairie
-- xmake require --info pour tchek les lib installed
-- xmake f -m debug or run pour changer de mode

