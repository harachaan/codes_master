### Folder の説明
修士の研究
- `estimation`
- `data`
- `gpryui`
- `hara_functions`

で完結してる．

以下はB4の研究なのでいらないっちゃいらない
- `regression`
- `create_traindata`

---
- estimation
    - `static_ver_surface`: 全部入ってる
    - `results`: 結果を出力するフォルダ

    ぶっちゃけほかはいらない．はず

    - `static_ver_surface`の中は，基本的にmainとgenがファイル名の最初にあるやつを動かす
        - `mainGPUKF_cmp_staticFixedPastMean.mlx`: 相対位置固定ver.
        - `mainGPUKF_cmp_staticPastMean_onOrbit.mlx`: 姿勢・軌道運動考慮ver.

    - codeの回す順番は上記の2つのファイルに書いてある．