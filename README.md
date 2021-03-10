# 概要
c++を使った電子状態計算ライブラリ

# 目的
- 計算手法の実装を考える/深い理解を得る
- c++に慣れる
- git/githubに慣れる
- API/Documentを書く

# To-do
- [x] xyzファイルから読み込み
- [x] 分子の構成
- [x] 基底関数の構成
- [x] Integral libraryの実装(libintを使って)
- [x] Restrected-Hartree-Fockのナイーブな実装
- [ ] Initial Density Matrixの構成
- [ ] DIISの実装
- [ ] UHFの実装
- [ ] Hatree-Fockの微分
- [ ] 核座標のOptimization
- [ ] (並列化?)
- [ ] CISDの実装
- [ ] [libxcを使ったDFT](https://www.tddft.org/programs/libxc/)

# Require
- [Libint(2.6.0)](https://github.com/evaleev/libint)
- [(Libxc)](https://www.tddft.org/programs/libxc/)