@echo off
setlocal enabledelayedexpansion
set mfile=%1

if "%mfile%"=="" (
    set mfile=Makefile
    call :main
) else if "%mfile%"=="-h" (
    call :help
) else call :main

:main
set /a v=0
set "all_exe="
for /f "tokens=*" %%i in ('type %mfile%') do (
    for %%j in (%%i) do (
        if %%j==install: (
            set /a v+=1
        )
        if !v!==1 (
            set "all_exe=!all_exe! %%j"
        )
        if %%j==all: (
            set /a v+=1
        )
    )
    if !v!==2 break
)
set "objects="
set exenum=0
for %%e in (%all_exe%) do (
    set v=0
    set /a exenum+=1
    for /f "tokens=*" %%i in ('type %mfile%') do (
        for %%j in (%%i) do (
            if !v!==1 (
                set "temp=%%j"
                set "objects=!objects! %%e !temp:~2,-1!"
                set /a v+=1
            )
            if %%j==%%e: (
                set /a v+=1
            )
        )
        if !v!==2 break
    )
)

echo Total %exenum% executable programs found in Makefile "%mfile%":%all_exe%
set v=0
set "exe="
for %%o in (%objects%) do (
    if !v!==0 (
        set "v=1"
        set "exe=%%o"

    ) else (
        echo Now working on !exe!
        set "v=0"
        set /a flag=0
        set "all_objects="
        for /f "tokens=*" %%i in ('type %mfile%') do (
            for %%j in (%%i) do (
                if !flag!==1 if %%j==include (
                    set /a flag+=1
                )
                if !flag!==1 (
                    set "all_objects=!all_objects! %%j"
                )
                if %%j==%%o (
                    set /a flag+=1
                )
            )
            if !flag!==2 break
        ) 
        echo    Writing !exe!.sln
        echo.> !exe!.sln
        echo Microsoft Visual Studio Solution File, Format Version 12.00>> !exe!.sln
        echo # Visual Studio 2013>> !exe!.sln
        echo VisualStudioVersion = 12.0.21005.1>> !exe!.sln
        echo MinimumVisualStudioVersion = 10.0.40219.1>> !exe!.sln
        echo Project^("{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}"^) = "!exe!", "!exe!.vcxproj", "{28E30FB4-9F17-4F20-985C-FF9E3F4E8B72}">> !exe!.sln
        echo EndProject>> !exe!.sln
        echo Global>> !exe!.sln
        echo 	GlobalSection^(SolutionConfigurationPlatforms^) = preSolution>> !exe!.sln
        echo 		Release^|x64 = Release^|x64>> !exe!.sln
        echo 	EndGlobalSection>> !exe!.sln
        echo 	GlobalSection^(ProjectConfigurationPlatforms^) = postSolution>> !exe!.sln
        echo 		{28E30FB4-9F17-4F20-985C-FF9E3F4E8B72}.Release^|x64.ActiveCfg = Release^|x64>> !exe!.sln
        echo 		{28E30FB4-9F17-4F20-985C-FF9E3F4E8B72}.Release^|x64.Build.0 = Release^|x64>> !exe!.sln
        echo 	EndGlobalSection>> !exe!.sln
        echo 	GlobalSection^(SolutionProperties^) = preSolution>> !exe!.sln
        echo 		HideSolutionNode = FALSE>> !exe!.sln
        echo 	EndGlobalSection>> !exe!.sln
        echo EndGlobal>> !exe!.sln
        
        echo    Writing !exe!.vcxproj
        echo ^<^?xml version="1.0" encoding="utf-8"^?^>> !exe!.vcxproj
        echo ^<Project DefaultTargets="Build" ToolsVersion="4.0" xmlns="http://schemas.microsoft.com/developer/msbuild/2003"^>>> !exe!.vcxproj
        echo   ^<ItemGroup Label="ProjectConfigurations"^>>> !exe!.vcxproj
        echo     ^<ProjectConfiguration Include="Release|x64"^>>> !exe!.vcxproj
        echo       ^<Configuration^>Release^</Configuration^>>> !exe!.vcxproj
        echo       ^<Platform^>x64^</Platform^>>> !exe!.vcxproj
        echo     ^</ProjectConfiguration^>>> !exe!.vcxproj
        echo   ^</ItemGroup^>>> !exe!.vcxproj
        echo   ^<PropertyGroup Label="Globals"^>>> !exe!.vcxproj
        echo     ^<ProjectGuid^>{28E30FB4-9F17-4F20-985C-FF9E3F4E8B72}^</ProjectGuid^>>> !exe!.vcxproj
        echo     ^<RootNamespace^>!exe!^</RootNamespace^>>> !exe!.vcxproj
        echo   ^</PropertyGroup^>>> !exe!.vcxproj
        echo   ^<Import Project="$(VCTargetsPath)\Microsoft.Cpp.Default.props" /^>>> !exe!.vcxproj
        echo   ^<PropertyGroup Condition="'$(Configuration)|$(Platform)'=='Release|x64'" Label="Configuration"^>>> !exe!.vcxproj
        echo     ^<ConfigurationType^>Application^</ConfigurationType^>>> !exe!.vcxproj
        echo     ^<UseDebugLibraries^>false^</UseDebugLibraries^>>> !exe!.vcxproj
        echo     ^<WholeProgramOptimization^>true^</WholeProgramOptimization^>>> !exe!.vcxproj
        echo     ^<CharacterSet^>MultiByte^</CharacterSet^>>> !exe!.vcxproj
        echo     ^<PlatformToolset^>v120^</PlatformToolset^>>> !exe!.vcxproj
        echo   ^</PropertyGroup^>>> !exe!.vcxproj
        echo   ^<Import Project="$(VCTargetsPath)\Microsoft.Cpp.props" /^>>> !exe!.vcxproj
        echo   ^<ImportGroup Label="ExtensionSettings"^>>> !exe!.vcxproj
        echo     ^<Import Project="$(VCTargetsPath)\BuildCustomizations\CUDA 10.2.props" /^>>> !exe!.vcxproj
        echo   ^</ImportGroup^>>> !exe!.vcxproj
        echo   ^<ImportGroup Label="PropertySheets" Condition="'$(Configuration)|$(Platform)'=='Release|x64'"^>>> !exe!.vcxproj
        echo     ^<Import Project="$(UserRootDir)\Microsoft.Cpp.$(Platform).user.props" Condition="exists('$(UserRootDir)\Microsoft.Cpp.$(Platform).user.props')" Label="LocalAppDataPlatform" /^>>> !exe!.vcxproj
        echo   ^</ImportGroup^>>> !exe!.vcxproj
        echo   ^<PropertyGroup Label="UserMacros" /^>>> !exe!.vcxproj
        echo   ^<ItemDefinitionGroup Condition="'$(Configuration)|$(Platform)'=='Release|x64'"^>>> !exe!.vcxproj
        echo     ^<ClCompile^>>> !exe!.vcxproj
        echo       ^<WarningLevel^>Level3^</WarningLevel^>>> !exe!.vcxproj
        echo       ^<Optimization^>MaxSpeed^</Optimization^>>> !exe!.vcxproj
        echo       ^<FunctionLevelLinking^>true^</FunctionLevelLinking^>>> !exe!.vcxproj
        echo       ^<IntrinsicFunctions^>true^</IntrinsicFunctions^>>> !exe!.vcxproj
        echo       ^<PreprocessorDefinitions^>WIN32;WIN64;NDEBUG;_CONSOLE;%%^(PreprocessorDefinitions^)^</PreprocessorDefinitions^>>> !exe!.vcxproj
        echo     ^</ClCompile^>>> !exe!.vcxproj
        echo     ^<Link^>>> !exe!.vcxproj
        echo       ^<GenerateDebugInformation^>true^</GenerateDebugInformation^>>> !exe!.vcxproj
        echo       ^<EnableCOMDATFolding^>true^</EnableCOMDATFolding^>>> !exe!.vcxproj
        echo       ^<OptimizeReferences^>true^</OptimizeReferences^>>> !exe!.vcxproj
        echo       ^<SubSystem^>Console^</SubSystem^>>> !exe!.vcxproj
        echo       ^<AdditionalDependencies^>cudadevrt.lib;cufft.lib;cudart_static.lib;kernel32.lib;user32.lib;gdi32.lib;winspool.lib;comdlg32.lib;advapi32.lib;shell32.lib;ole32.lib;oleaut32.lib;uuid.lib;odbc32.lib;odbccp32.lib;%%^(AdditionalDependencies^)^</AdditionalDependencies^>>> !exe!.vcxproj
        echo     ^</Link^>>> !exe!.vcxproj
        echo     ^<CudaCompile^>>> !exe!.vcxproj
        echo       ^<TargetMachinePlatform^>64^</TargetMachinePlatform^>>> !exe!.vcxproj
        echo       ^<GenerateRelocatableDeviceCode^>true^</GenerateRelocatableDeviceCode^>>> !exe!.vcxproj
        echo       ^<FastMath^>true^</FastMath^>>> !exe!.vcxproj
        echo       ^<Optimization^>O3^</Optimization^>>> !exe!.vcxproj
        echo       ^<CodeGeneration^>compute_50,sm_50^</CodeGeneration^>>> !exe!.vcxproj
        echo     ^</CudaCompile^>>> !exe!.vcxproj
        echo   ^</ItemDefinitionGroup^>>> !exe!.vcxproj
        echo   ^<ItemGroup^>>> !exe!.vcxproj
        for %%a in (!all_objects!) do (
            set temp=%%a
            set "temp=!temp:/=\!"
            echo     ^<CudaCompile Include="!temp:.o=.cu!" /^>>> !exe!.vcxproj
            echo     ^<ClInclude Include="!temp:.o=.cuh!" /^>>> !exe!.vcxproj
        )
        echo   ^</ItemGroup^>>> !exe!.vcxproj
        echo   ^<Import Project="$(VCTargetsPath)\Microsoft.Cpp.targets" /^>>> !exe!.vcxproj
        echo   ^<ImportGroup Label="ExtensionTargets"^>>> !exe!.vcxproj
        echo     ^<Import Project="$(VCTargetsPath)\BuildCustomizations\CUDA 10.2.targets" /^>>> !exe!.vcxproj
        echo   ^</ImportGroup^>>> !exe!.vcxproj
        echo ^</Project^>>> !exe!.vcxproj
        
        echo    Writing !exe!.vcxproj.filters
        echo ^<^?xml version="1.0" encoding="utf-8"^?^>> !exe!.vcxproj.filters
        echo ^<Project ToolsVersion="12.0" xmlns="http://schemas.microsoft.com/developer/msbuild/2003"^>>> !exe!.vcxproj.filters
        echo   ^<ItemGroup^>>> !exe!.vcxproj.filters
        for %%a in (!all_objects!) do (
            set temp=%%a
            set "temp=!temp:/=\!"
            set "ttemp=!temp:\=,!"
            set "filter="
            set tttemp=0
            for %%b in (!ttemp!) do (
                echo %%b| findstr \.o >nul || (
                    if !tttemp!==0 (set filter=%%b) else (set filter=!filter!\%%b)
                    set tttemp=1
                    echo     ^<Filter Include="!filter!" /^>>> !exe!.vcxproj.filters
                )                    
            )
            if !tttemp!==0 (
                echo     ^<CudaCompile Include="!temp:.o=.cu!" /^>>> !exe!.vcxproj.filters
                echo     ^<ClInclude Include="!temp:.o=.cuh!" /^>>> !exe!.vcxproj.filters
            ) else (
                echo     ^<CudaCompile Include="!temp:.o=.cu!"^>>> !exe!.vcxproj.filters
                echo       ^<Filter^>!filter!^</Filter^>>> !exe!.vcxproj.filters
                echo     ^</CudaCompile^>>> !exe!.vcxproj.filters
                echo     ^<ClInclude Include="!temp:.o=.cuh!"^>>> !exe!.vcxproj.filters
                echo       ^<Filter^>!filter!^</Filter^>>> !exe!.vcxproj.filters
                echo     ^</ClInclude^>>> !exe!.vcxproj.filters
            )
        )
        echo   ^</ItemGroup^>>> !exe!.vcxproj.filters
        echo ^</Project^>>> !exe!.vcxproj.filters
        
        echo    Writing !exe!.vcxproj.user
        echo ^<^?xml version="1.0" encoding="utf-8"^?^>> !exe!.vcxproj.user
        echo ^<Project ToolsVersion="12.0" xmlns="http://schemas.microsoft.com/developer/msbuild/2003"^>>> !exe!.vcxproj.user
        echo   ^<PropertyGroup Condition="'$(Configuration)|$(Platform)'=='Release|x64'"^>>> !exe!.vcxproj.user
        echo     ^<LocalDebuggerCommandArguments^>>> !exe!.vcxproj.user
        echo     ^</LocalDebuggerCommandArguments^>>> !exe!.vcxproj.user
        echo     ^<DebuggerFlavor^>WindowsLocalDebugger^</DebuggerFlavor^>>> !exe!.vcxproj.user
        echo     ^<LocalDebuggerWorkingDirectory^>$^(ProjectDir^)^</LocalDebuggerWorkingDirectory^>>> !exe!.vcxproj.user
        echo   ^</PropertyGroup^>>> !exe!.vcxproj.user
        echo   ^<PropertyGroup^>>> !exe!.vcxproj.user
        echo     ^<ShowAllFiles^>false^</ShowAllFiles^>>> !exe!.vcxproj.user
        echo   ^</PropertyGroup^>>> !exe!.vcxproj.user
        echo ^</Project^>>> !exe!.vcxproj.user
        
        
        echo.
    )
)
pause
exit

:help
echo .\vs_project_generator.bat [ARG]
echo    .\vs_project_generator.bat: generate visual studio project files from a Makefile
echo    ARG: the input parameter, it can be "-h", "--help", or the Makefile name. Default: Makefile
pause
exit