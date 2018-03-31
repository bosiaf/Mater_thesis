project(HIVEvo)

set(TARGET_NAME ${PROJECT_NAME})
set(headers
    epi.h
    global_fct.h
    par_class.h
    host.h
    strain.h
    )

    set(sources
        global_fct.cpp
        host.cpp
        par_class.cpp
        )
add_library(${TARGET_NAME} STATIC ${headers} ${sources})

target_include_directories(${TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

set_property(TARGET ${TARGET_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)

set_property(TARGET ${TARGET_NAME} PROPERTY FOLDER "HIVEvo")

hivevo_add_alias(${TARGET_NAME})

