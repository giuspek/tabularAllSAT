(set-logic QF_UF)
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
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(declare-fun A25 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A8))) (let ((.def_1 (or .def_0 A20))) (let ((.def_2 (not .def_1))) (let ((.def_3 (and A10 A5))) (let ((.def_4 (not .def_3))) (let ((.def_5 (and .def_4 .def_2))) (let ((.def_6 (not A26))) (let ((.def_7 (not A24))) (let ((.def_8 (or .def_7 .def_6))) (let ((.def_9 (not .def_8))) (let ((.def_10 (not A12))) (let ((.def_11 (not A19))) (let ((.def_12 (and .def_11 .def_10))) (let ((.def_13 (not .def_12))) (let ((.def_14 (or .def_13 .def_9))) (let ((.def_15 (not .def_14))) (let ((.def_16 (and .def_15 .def_5))) (let ((.def_17 (not .def_16))) (let ((.def_18 (not A3))) (let ((.def_19 (and A10 .def_18))) (let ((.def_20 (not A13))) (let ((.def_21 (not A23))) (let ((.def_22 (and .def_21 .def_20))) (let ((.def_23 (and .def_22 .def_19))) (let ((.def_24 (not .def_23))) (let ((.def_25 (not A11))) (let ((.def_26 (not A18))) (let ((.def_27 (= .def_26 .def_25))) (let ((.def_28 (not .def_27))) (let ((.def_29 (not A1))) (let ((.def_30 (and A20 .def_29))) (let ((.def_31 (not .def_30))) (let ((.def_32 (and .def_31 .def_28))) (let ((.def_33 (= .def_32 .def_24))) (let ((.def_34 (not .def_33))) (let ((.def_35 (or .def_34 .def_17))) (let ((.def_36 (or A5 A17))) (let ((.def_37 (not A7))) (let ((.def_38 (= .def_37 A19))) (let ((.def_39 (and .def_38 .def_36))) (let ((.def_40 (not .def_39))) (let ((.def_41 (not A14))) (let ((.def_42 (or .def_21 .def_41))) (let ((.def_43 (not .def_42))) (let ((.def_44 (and A7 A23))) (let ((.def_45 (not .def_44))) (let ((.def_46 (and .def_45 .def_43))) (let ((.def_47 (not .def_46))) (let ((.def_48 (or .def_47 .def_40))) (let ((.def_49 (not .def_48))) (let ((.def_50 (not A22))) (let ((.def_51 (or .def_21 .def_50))) (let ((.def_52 (or .def_37 .def_50))) (let ((.def_53 (= .def_52 .def_51))) (let ((.def_54 (not .def_53))) (let ((.def_55 (and .def_26 A19))) (let ((.def_56 (or A24 A20))) (let ((.def_57 (and .def_56 .def_55))) (let ((.def_58 (or .def_57 .def_54))) (let ((.def_59 (= .def_58 .def_49))) (let ((.def_60 (not .def_59))) (let ((.def_61 (= .def_60 .def_35))) (let ((.def_62 (not .def_61))) (let ((.def_63 (and A20 A28))) (let ((.def_64 (not .def_63))) (let ((.def_65 (and A11 A9))) (let ((.def_66 (not .def_65))) (let ((.def_67 (or .def_66 .def_64))) (let ((.def_68 (and A26 A8))) (let ((.def_69 (or A25 A22))) (let ((.def_70 (not .def_69))) (let ((.def_71 (or .def_70 .def_68))) (let ((.def_72 (not .def_71))) (let ((.def_73 (or .def_72 .def_67))) (let ((.def_74 (not .def_73))) (let ((.def_75 (not A27))) (let ((.def_76 (or .def_75 A25))) (let ((.def_77 (not .def_76))) (let ((.def_78 (not A28))) (let ((.def_79 (and .def_20 .def_78))) (let ((.def_80 (and .def_79 .def_77))) (let ((.def_81 (not .def_80))) (let ((.def_82 (or A25 .def_50))) (let ((.def_83 (or A14 .def_21))) (let ((.def_84 (and .def_83 .def_82))) (let ((.def_85 (not .def_84))) (let ((.def_86 (and .def_85 .def_81))) (let ((.def_87 (not .def_86))) (let ((.def_88 (and .def_87 .def_74))) (let ((.def_89 (and A18 .def_21))) (let ((.def_90 (not .def_89))) (let ((.def_91 (or .def_10 .def_7))) (let ((.def_92 (or .def_91 .def_90))) (let ((.def_93 (not .def_92))) (let ((.def_94 (not A29))) (let ((.def_95 (or A14 .def_94))) (let ((.def_96 (not .def_95))) (let ((.def_97 (not A25))) (let ((.def_98 (not A17))) (let ((.def_99 (or .def_98 .def_97))) (let ((.def_100 (not .def_99))) (let ((.def_101 (and .def_100 .def_96))) (let ((.def_102 (not .def_101))) (let ((.def_103 (= .def_102 .def_93))) (let ((.def_104 (not A9))) (let ((.def_105 (or .def_104 A14))) (let ((.def_106 (not .def_105))) (let ((.def_107 (or A3 A4))) (let ((.def_108 (and .def_107 .def_106))) (let ((.def_109 (not .def_108))) (let ((.def_110 (or A24 .def_0))) (let ((.def_111 (not .def_110))) (let ((.def_112 (or A6 .def_78))) (let ((.def_113 (not .def_112))) (let ((.def_114 (or .def_113 .def_111))) (let ((.def_115 (or .def_114 .def_109))) (let ((.def_116 (not .def_115))) (let ((.def_117 (= .def_116 .def_103))) (let ((.def_118 (not .def_117))) (let ((.def_119 (= .def_118 .def_88))) (let ((.def_120 (not .def_119))) (let ((.def_121 (or .def_120 .def_62))) .def_121)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)