set(dir ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND srcs "${dir}/api.f90" "${dir}/env.f90" "${dir}/main.f90")
set(srcs ${srcs} PARENT_SCOPE)
