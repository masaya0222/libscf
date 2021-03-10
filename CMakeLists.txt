cmake_minimum_required (VERSION 3.0)
project (LIBSCF CXX)

# executable のincludeディレクトリ
# set (CMAKE_POSITION_INDEPENDENT_CODE TRUE)
include_directories("src")
include_directories("/usr/local/libint/2.6.0/include")
include_directories("/usr/local/include/eigen3")
# Makeしていきたいディレクトリ
add_subdirectory("src")
# リンクしたいライブラリがある場所
link_libraries("/usr/local/libint/2.6.0/lib/libint2.a")

# 実行ファイルをつくる
add_executable(main main.cpp)
# 実行ファイルにshared libraryをリンク
target_link_libraries(main scf)