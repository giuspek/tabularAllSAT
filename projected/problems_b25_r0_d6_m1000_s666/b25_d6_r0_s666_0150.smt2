(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (not A0))) (let ((.def_1 (not A9))) (let ((.def_2 (or .def_1 .def_0))) (let ((.def_3 (not .def_2))) (let ((.def_4 (not A15))) (let ((.def_5 (not A18))) (let ((.def_6 (and .def_5 .def_4))) (let ((.def_7 (not .def_6))) (let ((.def_8 (or .def_7 .def_3))) (let ((.def_9 (not .def_8))) (let ((.def_10 (not A17))) (let ((.def_11 (or .def_1 .def_10))) (let ((.def_12 (not A20))) (let ((.def_13 (or A21 .def_12))) (let ((.def_14 (or .def_13 .def_11))) (let ((.def_15 (not .def_14))) (let ((.def_16 (and .def_15 .def_9))) (let ((.def_17 (not .def_16))) (let ((.def_18 (and A11 A18))) (let ((.def_19 (= A12 A15))) (let ((.def_20 (and .def_19 .def_18))) (let ((.def_21 (not A3))) (let ((.def_22 (or A0 .def_21))) (let ((.def_23 (or A16 A7))) (let ((.def_24 (not .def_23))) (let ((.def_25 (or .def_24 .def_22))) (let ((.def_26 (not .def_25))) (let ((.def_27 (and .def_26 .def_20))) (let ((.def_28 (not .def_27))) (let ((.def_29 (and .def_28 .def_17))) (let ((.def_30 (not .def_29))) (let ((.def_31 (and A6 A21))) (let ((.def_32 (not .def_31))) (let ((.def_33 (and A3 A12))) (let ((.def_34 (not .def_33))) (let ((.def_35 (and .def_34 .def_32))) (let ((.def_36 (not .def_35))) (let ((.def_37 (not A16))) (let ((.def_38 (and .def_37 A23))) (let ((.def_39 (not A8))) (let ((.def_40 (or .def_12 .def_39))) (let ((.def_41 (not .def_40))) (let ((.def_42 (and .def_41 .def_38))) (let ((.def_43 (and .def_42 .def_36))) (let ((.def_44 (and .def_37 A24))) (let ((.def_45 (not A22))) (let ((.def_46 (and A7 .def_45))) (let ((.def_47 (not .def_46))) (let ((.def_48 (or .def_47 .def_44))) (let ((.def_49 (not .def_48))) (let ((.def_50 (and A18 .def_10))) (let ((.def_51 (not .def_50))) (let ((.def_52 (or A11 A22))) (let ((.def_53 (not .def_52))) (let ((.def_54 (or .def_53 .def_51))) (let ((.def_55 (not .def_54))) (let ((.def_56 (and .def_55 .def_49))) (let ((.def_57 (not .def_56))) (let ((.def_58 (and .def_57 .def_43))) (let ((.def_59 (and .def_58 .def_30))) (let ((.def_60 (not .def_59))) (let ((.def_61 (not A4))) (let ((.def_62 (and A10 .def_61))) (let ((.def_63 (not .def_62))) (let ((.def_64 (or A22 A9))) (let ((.def_65 (or .def_64 .def_63))) (let ((.def_66 (not A23))) (let ((.def_67 (not A1))) (let ((.def_68 (and .def_67 .def_66))) (let ((.def_69 (or A10 .def_21))) (let ((.def_70 (or .def_69 .def_68))) (let ((.def_71 (not .def_70))) (let ((.def_72 (and .def_71 .def_65))) (let ((.def_73 (not .def_72))) (let ((.def_74 (or A5 .def_12))) (let ((.def_75 (not A5))) (let ((.def_76 (or .def_75 A6))) (let ((.def_77 (not .def_76))) (let ((.def_78 (or .def_77 .def_74))) (let ((.def_79 (not .def_78))) (let ((.def_80 (not A24))) (let ((.def_81 (and A5 .def_80))) (let ((.def_82 (not .def_81))) (let ((.def_83 (or .def_10 A23))) (let ((.def_84 (or .def_83 .def_82))) (let ((.def_85 (and .def_84 .def_79))) (let ((.def_86 (not .def_85))) (let ((.def_87 (and .def_86 .def_73))) (let ((.def_88 (not A11))) (let ((.def_89 (and .def_88 .def_37))) (let ((.def_90 (not .def_89))) (let ((.def_91 (or A9 .def_80))) (let ((.def_92 (not .def_91))) (let ((.def_93 (and .def_92 .def_90))) (let ((.def_94 (or A21 .def_67))) (let ((.def_95 (not .def_94))) (let ((.def_96 (not A7))) (let ((.def_97 (or A22 .def_96))) (let ((.def_98 (not .def_97))) (let ((.def_99 (and .def_98 .def_95))) (let ((.def_100 (and .def_99 .def_93))) (let ((.def_101 (or .def_12 A15))) (let ((.def_102 (and .def_45 .def_75))) (let ((.def_103 (not .def_102))) (let ((.def_104 (or .def_103 .def_101))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and .def_61 .def_67))) (let ((.def_107 (and A19 .def_66))) (let ((.def_108 (not .def_107))) (let ((.def_109 (= .def_108 .def_106))) (let ((.def_110 (not .def_109))) (let ((.def_111 (and .def_110 .def_105))) (let ((.def_112 (and .def_111 .def_100))) (let ((.def_113 (not .def_112))) (let ((.def_114 (and .def_113 .def_87))) (let ((.def_115 (or .def_114 .def_60))) (let ((.def_116 (not .def_115))) .def_116))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
