(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
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
(assert (let ((.def_0 (not A3))) (let ((.def_1 (not A16))) (let ((.def_2 (or .def_1 .def_0))) (let ((.def_3 (or A25 A23))) (let ((.def_4 (or .def_3 .def_2))) (let ((.def_5 (not .def_4))) (let ((.def_6 (not A14))) (let ((.def_7 (and A27 .def_6))) (let ((.def_8 (not .def_7))) (let ((.def_9 (or A29 A9))) (let ((.def_10 (not .def_9))) (let ((.def_11 (= .def_10 .def_8))) (let ((.def_12 (and .def_11 .def_5))) (let ((.def_13 (not A20))) (let ((.def_14 (not A5))) (let ((.def_15 (and .def_14 .def_13))) (let ((.def_16 (not .def_15))) (let ((.def_17 (or .def_14 A25))) (let ((.def_18 (not .def_17))) (let ((.def_19 (or .def_18 .def_16))) (let ((.def_20 (not .def_19))) (let ((.def_21 (not A7))) (let ((.def_22 (not A23))) (let ((.def_23 (and .def_22 .def_21))) (let ((.def_24 (not .def_23))) (let ((.def_25 (not A9))) (let ((.def_26 (and .def_25 A9))) (let ((.def_27 (not .def_26))) (let ((.def_28 (and .def_27 .def_24))) (let ((.def_29 (and .def_28 .def_20))) (let ((.def_30 (not .def_29))) (let ((.def_31 (and .def_30 .def_12))) (let ((.def_32 (not .def_31))) (let ((.def_33 (not A19))) (let ((.def_34 (not A6))) (let ((.def_35 (= .def_34 .def_33))) (let ((.def_36 (not .def_35))) (let ((.def_37 (not A17))) (let ((.def_38 (or .def_37 A10))) (let ((.def_39 (and .def_38 .def_36))) (let ((.def_40 (not .def_39))) (let ((.def_41 (not A26))) (let ((.def_42 (and .def_41 A23))) (let ((.def_43 (not .def_42))) (let ((.def_44 (not A24))) (let ((.def_45 (= A21 .def_44))) (let ((.def_46 (not .def_45))) (let ((.def_47 (or .def_46 .def_43))) (let ((.def_48 (not .def_47))) (let ((.def_49 (or .def_48 .def_40))) (let ((.def_50 (not A21))) (let ((.def_51 (and .def_50 A25))) (let ((.def_52 (not .def_51))) (let ((.def_53 (not A27))) (let ((.def_54 (or .def_1 .def_53))) (let ((.def_55 (or .def_54 .def_52))) (let ((.def_56 (not .def_55))) (let ((.def_57 (not A15))) (let ((.def_58 (or A19 .def_57))) (let ((.def_59 (not .def_58))) (let ((.def_60 (and A13 .def_57))) (let ((.def_61 (not .def_60))) (let ((.def_62 (and .def_61 .def_59))) (let ((.def_63 (not .def_62))) (let ((.def_64 (or .def_63 .def_56))) (let ((.def_65 (not .def_64))) (let ((.def_66 (or .def_65 .def_49))) (let ((.def_67 (or .def_66 .def_32))) (let ((.def_68 (not A12))) (let ((.def_69 (not A28))) (let ((.def_70 (and .def_69 .def_68))) (let ((.def_71 (not A29))) (let ((.def_72 (not A10))) (let ((.def_73 (and .def_72 .def_71))) (let ((.def_74 (or .def_73 .def_70))) (let ((.def_75 (or .def_13 .def_34))) (let ((.def_76 (and A19 A18))) (let ((.def_77 (= .def_76 .def_75))) (let ((.def_78 (not .def_77))) (let ((.def_79 (or .def_78 .def_74))) (let ((.def_80 (not A13))) (let ((.def_81 (= .def_13 .def_80))) (let ((.def_82 (not .def_81))) (let ((.def_83 (and A4 A16))) (let ((.def_84 (and .def_83 .def_82))) (let ((.def_85 (not .def_84))) (let ((.def_86 (or A15 .def_80))) (let ((.def_87 (not .def_86))) (let ((.def_88 (not A22))) (let ((.def_89 (and A8 .def_88))) (let ((.def_90 (not .def_89))) (let ((.def_91 (and .def_90 .def_87))) (let ((.def_92 (not .def_91))) (let ((.def_93 (and .def_92 .def_85))) (let ((.def_94 (not .def_93))) (let ((.def_95 (or .def_94 .def_79))) (let ((.def_96 (not .def_95))) (let ((.def_97 (and .def_13 .def_34))) (let ((.def_98 (or A13 A9))) (let ((.def_99 (or .def_98 .def_97))) (let ((.def_100 (not .def_99))) (let ((.def_101 (or A15 A28))) (let ((.def_102 (not .def_101))) (let ((.def_103 (and A0 A4))) (let ((.def_104 (or .def_103 .def_102))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and .def_105 .def_100))) (let ((.def_107 (and .def_53 A12))) (let ((.def_108 (not A8))) (let ((.def_109 (= .def_88 .def_108))) (let ((.def_110 (or .def_109 .def_107))) (let ((.def_111 (not .def_110))) (let ((.def_112 (not A4))) (let ((.def_113 (not A25))) (let ((.def_114 (or .def_113 .def_112))) (let ((.def_115 (not .def_114))) (let ((.def_116 (and .def_112 .def_14))) (let ((.def_117 (and .def_116 .def_115))) (let ((.def_118 (not .def_117))) (let ((.def_119 (or .def_118 .def_111))) (let ((.def_120 (not .def_119))) (let ((.def_121 (and .def_120 .def_106))) (let ((.def_122 (not .def_121))) (let ((.def_123 (or .def_122 .def_96))) (let ((.def_124 (not .def_123))) (let ((.def_125 (and .def_124 .def_67))) (let ((.def_126 (not .def_125))) .def_126))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
