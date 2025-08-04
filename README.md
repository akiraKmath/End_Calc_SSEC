<!-- (仮)日本語版  
# 超特異楕円曲線の自己準同型環計算アルゴリズム
有限体 $\mathbb{F}_q$上の超特異楕円曲線 $E$の自己準同型環 $\mathrm{End}(E)$の基底計算のコードを紹介する.  
本コードは~~で提案されたものの実装したものである.
## 理論的説明
詳細は~~を参照.
## 実装の説明
本コードは以下の3つのステップに分けることができる:
1. End環計算:
    1. 同種写像サイクル探索:
    1. 四元数代数の元への変換:
    1. 極大整環であるかの判定:
2. KLPTアルゴリズム計算:
3. Deuring対応の計算

## テストコードについて
テストコードは`sage test_end_calc.sage`で実行可能である.   
編集可能なパラメータは以下の通りである:
- `k`: $\mathbb{F}_q$の標数 $p$のビット数
- `ells_num`: サイクルの次数に用いる素数の数
- `collect_nums`: 2,3-isogenyの数
- `Fp_defined`: 基礎体(`Fp=1`なら素体, `Fp=2`なら拡大体, `Fp=0`なら制限なし)
- `D`: 
- `is_elkies`: Elkies素数の利用の有無
- `test_times`: 実行回数 -->

# Algorithm for computing the endomorphism ring of a supersingular elliptic curve over a finite field
A SageMath implementation of the algorithm following the paper "Implementation report on computing the endomorphism ring of a supersingular elliptic curve, by Yuta Kambe, Akira Katayama, Kazuki Komine, Yusuke Aikawa, Yuki Ishihara, Masaya Yasuda, and Kazuhiro Yokoyama.
## Theoretical Background

## The details of the implementation
The alogrithm can be divided into the following steps:
1. Computing the endomorphism ring:
    1. Finding an isogeny cycle on $E$
    1. Conversion to an element of $\mathbb{Q}$-quanternion algebra
    1. Constructing a basis of the endomorphism ring as a lattice
2. Computing the KLPT algorithm
3. Computing the Deuring correspondence

## Example Usage
This implementation requires the database of modular polynomials.
Download the database from [https://aur.archlinux.org/packages/sage-data-kohel](https://aur.archlinux.org/packages/sage-data-kohel) before running the code.
You can run `sage example_code.sage`
The following parameters in `example_code.sage` are editable:
- `k`: The bit number of $p$ the characteristic of $\mathbb{F}_q$
- `ells_num`: The number of prime numbers for finding isogeny cycles 
- `collect_nums`: The number of collected 2,3-isogenies
- `Fp_defined`: The base field of $E$(if `Fp=1`, prime field, if `Fp=2`, the extension field (the default), if `Fp=0`, either of the prime field or the extension field)
- `is_elkies`: Either using Elkies' method or not
- `test_times`: The number of executions

## License
An implementation program of computing the endomorphism ring of a given supersingular elliptic curve besed on finding isogeny cycles.
(C) 2025 Mitsubisi Electric, Rikkyo University, Created by Yuta Kambe, Akira Katayama, Kazuki Komine, Yusuke Aikawa, Yuki Ishihara, Masaya Yasuda, Kazuhiro Yokoyama.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
