(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
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
(declare-fun A25 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (and A5 A19))) (let ((.def_1 (or A20 A10))) (let ((.def_2 (not .def_1))) (let ((.def_3 (and .def_2 .def_0))) (let ((.def_4 (or A8 A4))) (let ((.def_5 (not A25))) (let ((.def_6 (not A4))) (let ((.def_7 (or .def_6 .def_5))) (let ((.def_8 (not .def_7))) (let ((.def_9 (and .def_8 .def_4))) (let ((.def_10 (and .def_9 .def_3))) (let ((.def_11 (not A22))) (let ((.def_12 (not A3))) (let ((.def_13 (or .def_12 .def_11))) (let ((.def_14 (not .def_13))) (let ((.def_15 (not A1))) (let ((.def_16 (or A0 .def_15))) (let ((.def_17 (not .def_16))) (let ((.def_18 (or .def_17 .def_14))) (let ((.def_19 (not .def_18))) (let ((.def_20 (not A12))) (let ((.def_21 (= .def_20 A23))) (let ((.def_22 (not .def_21))) (let ((.def_23 (not A24))) (let ((.def_24 (not A28))) (let ((.def_25 (= .def_24 .def_23))) (let ((.def_26 (not .def_25))) (let ((.def_27 (or .def_26 .def_22))) (let ((.def_28 (= .def_27 .def_19))) (let ((.def_29 (not .def_28))) (let ((.def_30 (= .def_29 .def_10))) (let ((.def_31 (not .def_30))) (let ((.def_32 (not A9))) (let ((.def_33 (= .def_32 A7))) (let ((.def_34 (not .def_33))) (let ((.def_35 (or A13 A4))) (let ((.def_36 (not .def_35))) (let ((.def_37 (= .def_36 .def_34))) (let ((.def_38 (not .def_37))) (let ((.def_39 (not A29))) (let ((.def_40 (not A17))) (let ((.def_41 (or .def_40 .def_39))) (let ((.def_42 (not .def_41))) (let ((.def_43 (or A7 A10))) (let ((.def_44 (not .def_43))) (let ((.def_45 (and .def_44 .def_42))) (let ((.def_46 (not .def_45))) (let ((.def_47 (or .def_46 .def_38))) (let ((.def_48 (and .def_39 A23))) (let ((.def_49 (not .def_48))) (let ((.def_50 (not A15))) (let ((.def_51 (or A23 .def_50))) (let ((.def_52 (not .def_51))) (let ((.def_53 (and .def_52 .def_49))) (let ((.def_54 (not A20))) (let ((.def_55 (and .def_20 .def_54))) (let ((.def_56 (not .def_55))) (let ((.def_57 (not A7))) (let ((.def_58 (and A8 .def_57))) (let ((.def_59 (not .def_58))) (let ((.def_60 (or .def_59 .def_56))) (let ((.def_61 (or .def_60 .def_53))) (let ((.def_62 (not .def_61))) (let ((.def_63 (= .def_62 .def_47))) (let ((.def_64 (or .def_63 .def_31))) (let ((.def_65 (and A16 A17))) (let ((.def_66 (not A18))) (let ((.def_67 (= .def_32 .def_66))) (let ((.def_68 (or .def_67 .def_65))) (let ((.def_69 (not .def_68))) (let ((.def_70 (not A27))) (let ((.def_71 (and .def_70 .def_40))) (let ((.def_72 (not .def_71))) (let ((.def_73 (= A24 A8))) (let ((.def_74 (not .def_73))) (let ((.def_75 (and .def_74 .def_72))) (let ((.def_76 (not .def_75))) (let ((.def_77 (and .def_76 .def_69))) (let ((.def_78 (not .def_77))) (let ((.def_79 (or A15 .def_5))) (let ((.def_80 (or A4 .def_20))) (let ((.def_81 (or .def_80 .def_79))) (let ((.def_82 (not A21))) (let ((.def_83 (not A11))) (let ((.def_84 (and .def_83 .def_82))) (let ((.def_85 (not .def_84))) (let ((.def_86 (not A26))) (let ((.def_87 (or A5 .def_86))) (let ((.def_88 (and .def_87 .def_85))) (let ((.def_89 (not .def_88))) (let ((.def_90 (or .def_89 .def_81))) (let ((.def_91 (not .def_90))) (let ((.def_92 (or .def_91 .def_78))) (let ((.def_93 (not .def_92))) (let ((.def_94 (not A10))) (let ((.def_95 (and A9 .def_94))) (let ((.def_96 (not .def_95))) (let ((.def_97 (and .def_39 .def_50))) (let ((.def_98 (not .def_97))) (let ((.def_99 (or .def_98 .def_96))) (let ((.def_100 (or .def_12 A18))) (let ((.def_101 (not .def_100))) (let ((.def_102 (and .def_11 A26))) (let ((.def_103 (and .def_102 .def_101))) (let ((.def_104 (not .def_103))) (let ((.def_105 (and .def_104 .def_99))) (let ((.def_106 (and A7 A0))) (let ((.def_107 (not .def_106))) (let ((.def_108 (not A5))) (let ((.def_109 (or .def_108 A9))) (let ((.def_110 (not .def_109))) (let ((.def_111 (and .def_110 .def_107))) (let ((.def_112 (not .def_111))) (let ((.def_113 (not A14))) (let ((.def_114 (= .def_113 A14))) (let ((.def_115 (or A14 .def_24))) (let ((.def_116 (not .def_115))) (let ((.def_117 (or .def_116 .def_114))) (let ((.def_118 (= .def_117 .def_112))) (let ((.def_119 (or .def_118 .def_105))) (let ((.def_120 (and .def_119 .def_93))) (let ((.def_121 (not .def_120))) (let ((.def_122 (and .def_121 .def_64))) .def_122))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
