# builder

## 圧力をかけ始めたらエラー
- 状況：charmmguiの拘束力でstep6.4(1fs→2fs)で平衡化してたらエラー
- 結論：charmmguiほどpackingがうまくないので、tleapで`setbox SYS vdw`とする必要がある+平衡化時の拘束力を弱める必要がある

## acpype grompp segmentation fault without restraints
- 状況：acpypeで作った系に対して拘束力を弱めてもproductionMDのgromppでセグフォ
- 結論：gromppのbugぽい、もしくは非推奨
- 回避策：生産MDでも拘束力=0でMDする
- 作業メモ：
  - ifdefをtopの中に入れたらどうなる？
    - 同じエラー
  - itpをまとめて一緒にしたらどうなる？
    - 同じエラー
  - 脂質への拘束を無くしたらどうなる？
    - 脂質への拘束を無くしたらちゃんとずれてるので多分拘束は効いてる
    - 脂質への拘束無くしたら動いた。。。(拘束が一つだけだとOK)
  - parmedで作った系に位置拘束を2個付いれたらどうなる？
    - charmmguiのやつは2つ入れても動く。。。、多分moleculetypeを分けているから
  - acpypeで他の系を作ったらどうなる？
    - リガンドを揃えないと読み込んで貰えない
    - ダメだった、、

