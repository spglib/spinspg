from spinspg.pointgroup import SPIN_POINT_GROUP_TYPES


def main():
    symbols = {
        "1": "1",
        "-1": r"\overline{1}",
        "2": "2",
        "m": "m",
        "2/m": "2/m",
        "222": "222",
        "mm2": "mm2",
        "mmm": "mmm",
        "4": "4",
        "-4": r"\overline{4}",
        "4/m": "4/m",
        "422": "422",
        "4mm": "4mm",
        "-42m": r"\overline{4}2m",
        "4/mmm": "4/mmm",
        "3": "3",
        "-3": r"\overline{3}",
        "32": "32",
        "3m": "3m",
        "-3m": r"\overline{3}m",
        "6": "6",
        "-6": r"\overline{6}",
        "6/m": "6/m",
        "622": "622",
        "6mm": "6mm",
        "-6m2": r"\overline{6}m2",
        "6/mmm": "6/mmm",
        "23": "23",
        "m-3": r"m\overline{3}",
        "432": "432",
        "-43m": r"\overline{4}3m",
        "m-3m": r"m\overline{3}m",
    }
    generator_simple_labels = {
        "1": ["1"],
        "-1": [r"\overline{1}"],
        "2": ["2"],
        "m": ["m"],
        "2/m": ["2", "m"],
        "222": ["2", "2", "2"],
        "mm2": ["m", "m", "2"],
        "mmm": ["m", "m", "m"],
        "4": ["4"],
        "-4": [r"\overline{4}"],
        "4/m": ["4", "m"],
        "422": ["4", "2", "2"],
        "4mm": ["4", "m", "m"],
        "-42m": [r"\overline{4}", "2", "m"],
        "4/mmm": ["4", "m", "m", "m"],
        "3": ["3"],
        "-3": [r"\overline{3}"],
        "32": ["3", "2"],
        "3m": ["3", "m"],
        "-3m": [r"\overline{3}", "m"],
        "6": ["6"],
        "-6": [r"\overline{6}"],
        "6/m": ["6", "m"],
        "622": ["6", "2", "2"],
        "6mm": ["6", "m", "m"],
        "-6m2": [r"\overline{6}", "m", "2"],
        "6/mmm": ["6", "m", "m", "m"],
        "23": ["2", "3"],
        "m-3": ["m", r"\overline{3}"],
        "432": ["4", "3", "2"],
        "-43m": [r"\overline{4}", "3", "m"],
        "m-3m": ["m", r"\overline{3}", "m"],
    }
    operation_labels = {
        "1": ["1"],
        "-1": ["1", r"\overline{1}"],
        "2": ["1", "2"],
        "m": ["1", "m"],
        "2/m": ["1", "2", r"\overline{1}", "m"],
        "222": ["1", r"2_{001}", r"2_{010}", r"2_{100}"],
        "mm2": ["1", "2", r"m_{010}", r"m_{100}"],
        "mmm": [
            "1",
            r"2_{001}",
            r"2_{010}",
            r"2_{100}",
            r"\overline{1}",
            r"m_{001}",
            r"m_{010}",
            r"m_{100}",
        ],
        "4": ["1", "2", r"4^{+}", r"4^{-}"],
        "-4": ["1", "2", r"\overline{4}^{+}", r"\overline{4}^{-}"],
        "4/m": [
            "1",
            "2",
            r"4^{+}",
            r"4^{-}",
            r"\overline{1}",
            "m",
            r"\overline{4}^{+}",
            r"\overline{4}^{-}",
        ],
        "422": [
            "1",
            r"2_{001}",
            r"4^{+}",
            r"4^{-}",
            r"2_{010}",
            r"2_{100}",
            r"2_{110}",
            r"2_{1\overline{1}0}",
        ],
        "4mm": [
            "1",
            "2",
            r"4^{+}",
            r"4^{-}",
            r"m_{010}",
            r"m_{100}",
            r"m_{110}",
            r"m_{1\overline{1}0}",
        ],
        "-42m": [
            "1",
            r"2_{001}",
            r"\overline{4}^{+}",
            r"\overline{4}^{-}",
            r"2_{010}",
            r"2_{100}",
            r"m_{110}",
            r"m_{1\overline{1}0}",
        ],
        "4/mmm": [
            "1",
            r"2_{001}",
            r"4^{+}",
            r"4^{-}",
            r"2_{010}",
            r"2_{100}",
            r"2_{110}",
            r"2_{1\overline{1}0}",
            r"\overline{1}",
            r"m_{001}",
            r"\overline{4}^{+}",
            r"\overline{4}^{-}",
            r"m_{010}",
            r"m_{100}",
            r"m_{110}",
            r"m_{1\overline{1}0}",
        ],
        "3": [
            "1",
            r"3^{+}",
            r"3^{-}",
        ],
        "-3": [
            "1",
            r"3^{+}",
            r"3^{-}",
            r"\overline{1}",
            r"\overline{3}^{+}",
            r"\overline{3}^{-}",
        ],
        # 312
        "32": [
            "1",
            r"3^{+}",
            r"3^{-}",
            r"2_{1\overline{1}0}",
            r"2_{120}",
            r"2_{210}",
        ],
        # 3m1
        "3m": [
            "1",
            r"3^{+}",
            r"3^{-}",
            r"m_{110}",
            r"m_{100}",
            r"m_{010}",
        ],
        # -31m
        "-3m": [
            "1",
            r"3^{+}",
            r"3^{-}",
            r"2_{1\overline{1}0}",
            r"2_{120}",
            r"2_{210}",
            r"\overline{1}",
            r"\overline{3}^{+}",
            r"\overline{3}^{-}",
            r"m_{1\overline{1}0}",
            r"m_{120}",
            r"m_{210}",
        ],
        "6": [
            "1",
            r"3^{+}",
            r"3^{-}",
            "2",
            r"6^{-}",
            r"6^{+}",
        ],
        "-6": [
            "1",
            r"3^{+}",
            r"3^{-}",
            "m",
            r"\overline{6}^{-}",
            r"\overline{6}^{+}",
        ],
        "6/m": [
            "1",
            r"3^{+}",
            r"3^{-}",
            "2",
            r"6^{-}",
            r"6^{+}",
            r"\overline{1}",
            r"\overline{3}^{+}",
            r"\overline{3}^{-}",
            "m",
            r"\overline{6}^{-}",
            r"\overline{6}^{+}",
        ],
        "622": [
            "1",
            r"3^{+}",
            r"3^{-}",
            r"2_{001}",
            r"6^{-}",
            r"6^{+}",
            r"2_{110}",
            r"2_{100}",
            r"2_{010}",
            r"2_{1\overline{1}0}",
            r"2_{120}",
            r"2_{210}",
        ],
        "6mm": [
            "1",
            r"3^{+}",
            r"3^{-}",
            r"2_{001}",
            r"6^{-}",
            r"6^{+}",
            r"m_{110}",
            r"m_{100}",
            r"m_{010}",
            r"m_{1\overline{1}0}",
            r"m_{120}",
            r"m_{210}",
        ],
        "-6m2": [
            "1",
            r"3^{+}",
            r"3^{-}",
            r"m_{001}",
            r"\overline{6}^{-}",
            r"\overline{6}^{+}",
            r"m_{110}",
            r"m_{100}",
            r"m_{010}",
            r"2_{1\overline{1}0}",
            r"2_{120}",
            r"2_{210}",
        ],
        "6/mmm": [
            "1",
            r"3^{+}",
            r"3^{-}",
            r"2_{001}",
            r"6^{-}",
            r"6^{+}",
            r"2_{110}",
            r"2_{100}",
            r"2_{010}",
            r"2_{1\overline{1}0}",
            r"2_{120}",
            r"2_{210}",
            r"\overline{1}",
            r"\overline{3}^{+}",
            r"\overline{3}^{-}",
            r"m_{001}",
            r"\overline{6}^{-}",
            r"\overline{6}^{+}",
            r"m_{110}",
            r"m_{100}",
            r"m_{010}",
            r"m_{1\overline{1}0}",
            r"m_{120}",
            r"m_{210}",
        ],
        "23": [
            "1",
            r"2_{001}",
            r"2_{010}",
            r"2_{100}",
            r"3^{+}_{111}",
            r"3^{+}_{\overline{1}1\overline{1}}",
            r"3^{+}_{1\overline{1}\overline{1}}",
            r"3^{+}_{\overline{1}\overline{1}1}",
            r"3^{-}_{111}",
            r"3^{-}_{1\overline{1}\overline{1}}",
            r"3^{-}_{\overline{1}\overline{1}1}",
            r"3^{-}_{\overline{1}1\overline{1}}",
        ],
        "m-3": [
            "1",
            r"2_{001}",
            r"2_{010}",
            r"2_{100}",
            r"3^{+}_{111}",
            r"3^{+}_{\overline{1}1\overline{1}}",
            r"3^{+}_{1\overline{1}\overline{1}}",
            r"3^{+}_{\overline{1}\overline{1}1}",
            r"3^{-}_{111}",
            r"3^{-}_{1\overline{1}\overline{1}}",
            r"3^{-}_{\overline{1}\overline{1}1}",
            r"3^{-}_{\overline{1}1\overline{1}}",
            r"\overline{1}",
            r"m_{001}",
            r"m_{010}",
            r"m_{100}",
            r"\overline{3}^{+}_{111}",
            r"\overline{3}^{+}_{\overline{1}1\overline{1}}",
            r"\overline{3}^{+}_{1\overline{1}\overline{1}}",
            r"\overline{3}^{+}_{\overline{1}\overline{1}1}",
            r"\overline{3}^{-}_{111}",
            r"\overline{3}^{-}_{1\overline{1}\overline{1}}",
            r"\overline{3}^{-}_{\overline{1}\overline{1}1}",
            r"\overline{3}^{-}_{\overline{1}1\overline{1}}",
        ],
        "432": [
            "1",
            r"2_{001}",
            r"2_{010}",
            r"2_{100}",
            r"3^{+}_{111}",
            r"3^{+}_{\overline{1}1\overline{1}}",
            r"3^{+}_{1\overline{1}\overline{1}}",
            r"3^{+}_{\overline{1}\overline{1}1}",
            r"3^{-}_{111}",
            r"3^{-}_{1\overline{1}\overline{1}}",
            r"3^{-}_{\overline{1}\overline{1}1}",
            r"3^{-}_{\overline{1}1\overline{1}}",
            r"2_{110}",
            r"2_{1\overline{1}0}",
            r"4^{-}_{001}",
            r"4^{+}_{001}",
            r"4^{-}_{100}",
            r"2_{011}",
            r"2_{01\overline{1}}",
            r"4^{+}_{100}",
            r"4^{+}_{010}",
            r"2_{101}",
            r"4^{-}_{010}",
            r"2_{\overline{1}01}",
        ],
        "-43m": [
            "1",
            r"2_{001}",
            r"2_{010}",
            r"2_{100}",
            r"3^{+}_{111}",
            r"3^{+}_{\overline{1}1\overline{1}}",
            r"3^{+}_{1\overline{1}\overline{1}}",
            r"3^{+}_{\overline{1}\overline{1}1}",
            r"3^{-}_{111}",
            r"3^{-}_{1\overline{1}\overline{1}}",
            r"3^{-}_{\overline{1}\overline{1}1}",
            r"3^{-}_{\overline{1}1\overline{1}}",
            r"m_{1\overline{1}0}",
            r"m_{110}",
            r"\overline{4}^{-}_{001}",
            r"\overline{4}^{+}_{001}",
            r"m_{01\overline{1}}",
            r"\overline{4}^{-}_{100}",
            r"\overline{4}^{-}_{100}",
            r"m_{01\overline{1}}",
            r"m_{\overline{1}01}",
            r"\overline{4}^{-}_{010}",
            r"m_{101}",
            r"\overline{4}^{+}_{010}",
        ],
        "m-3m": [
            "1",
            r"2_{001}",
            r"2_{010}",
            r"2_{100}",
            r"3^{+}_{111}",
            r"3^{+}_{\overline{1}1\overline{1}}",
            r"3^{+}_{1\overline{1}\overline{1}}",
            r"3^{+}_{\overline{1}\overline{1}1}",
            r"3^{-}_{111}",
            r"3^{-}_{1\overline{1}\overline{1}}",
            r"3^{-}_{\overline{1}\overline{1}1}",
            r"3^{-}_{\overline{1}1\overline{1}}",
            r"2_{110}",
            r"2_{1\overline{1}0}",
            r"4^{-}_{001}",
            r"4^{+}_{001}",
            r"4^{-}_{100}",
            r"2_{011}",
            r"2_{01\overline{1}}",
            r"4^{+}_{100}",
            r"4^{+}_{010}",
            r"2_{101}",
            r"4^{-}_{010}",
            r"2_{\overline{1}01}",
            # 24-
            r"\overline{1}",
            r"m_{001}",
            r"m_{010}",
            r"m_{100}",
            r"\overline{3}^{+}_{111}",
            r"\overline{3}^{+}_{\overline{1}1\overline{1}}",
            r"\overline{3}^{+}_{1\overline{1}\overline{1}}",
            r"\overline{3}^{+}_{\overline{1}\overline{1}1}",
            r"\overline{3}^{-}_{111}",
            r"\overline{3}^{-}_{1\overline{1}\overline{1}}",
            r"\overline{3}^{-}_{\overline{1}\overline{1}1}",
            r"\overline{3}^{-}_{\overline{1}1\overline{1}}",
            r"m_{110}",
            r"m_{1\overline{1}0}",
            r"\overline{4}^{+}_{001}",
            r"\overline{4}^{+}_{001}" r"\overline{4}^{-}_{100}",
            r"m_{011}",
            r"m_{01\overline{1}}",
            r"\overline{4}^{+}_{100}",
            r"\overline{4}^{+}_{010}",
            r"m_{101}",
            r"\overline{4}^{-}_{010}",
            r"m_{\overline{1}01}",
        ],
    }

    # Assert operation_labels
    for pg_symbol, labels in operation_labels.items():
        if len(set(labels)) != len(labels):
            print(f"Wrong labels in {pg_symbol}")

    contents = []

    contents.append(
        r"""\section{\label{appx:spin_point_group_table}Tabulation of nontrivial spin point group types}

\newpage

\begin{longtable}{ccccc}
  \caption{Nontrivial spin point group types}
  \label{tab:spin_point_group_types} \\
  \hline \hline
  Litvin number & $R$ & $r$ & $B$ & Symbol \\
  \hline
  \endfirsthead
  \multicolumn{5}{c}{\tablename\ \thetable\ (\textit{cont.})} \\
  Litvin number & $R$ & $r$ & $B$ & Symbol \\
  \hline
  \endhead
  \endfoot
  \endlastfoot"""
    )

    for symbol_R, datum_R in SPIN_POINT_GROUP_TYPES.items():
        new_R = True
        for symbol_r, datum_R_r in datum_R.items():  # type: ignore
            new_r = True
            for symbol_B, datum_R_r_B in datum_R_r.items():
                new_B = True
                for number, mapping in datum_R_r_B:
                    print(f"#{number}: R={symbol_R}, r={symbol_r}, B={symbol_B}")

                    entries = [
                        "  " + str(number),
                    ]

                    if new_B and new_r and new_R:
                        entries.append("$" + symbols[symbol_r] + "$")
                    else:
                        entries.append("")

                    if new_B and new_r:
                        entries.append("$" + symbols[symbol_r] + "$")
                    else:
                        entries.append("")

                    if new_B:
                        entries.append("$" + symbols[symbol_B] + "$")
                    else:
                        entries.append("")

                    spin_operations = [operation_labels[symbol_B][idx] for idx in mapping]
                    spatial_operations = generator_simple_labels[symbol_R]
                    has_slash = "/" in symbol_R
                    entry = r"$"
                    for i, (rot, srot) in enumerate(zip(spatial_operations, spin_operations)):
                        if has_slash and (i == 1):
                            entry += r" / "
                        entry += r"{}^{" + srot + r"} " + rot + " "
                    entry += r"$"

                    entries.append(entry)
                    contents.append(r" & ".join(entries) + r"\\")

                    new_B = False
                new_r = False
            new_R = False

    contents.append(r"  \hline\hline")
    contents.append(r"\end{longtable}")

    with open("spin_point_group_table.tex", "w") as f:
        for line in contents:
            f.write(line)
            f.write("\n")


if __name__ == "__main__":
    main()
