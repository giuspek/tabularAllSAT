(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
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
(assert (let ((.def_0 (not A2))) (let ((.def_1 (and .def_0 A8))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A4))) (let ((.def_4 (and .def_3 A15))) (let ((.def_5 (and .def_4 .def_2))) (let ((.def_6 (not .def_5))) (let ((.def_7 (or A23 A15))) (let ((.def_8 (not .def_7))) (let ((.def_9 (and A5 A16))) (let ((.def_10 (not .def_9))) (let ((.def_11 (= .def_10 .def_8))) (let ((.def_12 (not .def_11))) (let ((.def_13 (or .def_12 .def_6))) (let ((.def_14 (not .def_13))) (let ((.def_15 (not A17))) (let ((.def_16 (or A15 .def_15))) (let ((.def_17 (not .def_16))) (let ((.def_18 (not A7))) (let ((.def_19 (not A23))) (let ((.def_20 (and .def_19 .def_18))) (let ((.def_21 (and .def_20 .def_17))) (let ((.def_22 (not A15))) (let ((.def_23 (and A1 .def_22))) (let ((.def_24 (not A20))) (let ((.def_25 (not A12))) (let ((.def_26 (or .def_25 .def_24))) (let ((.def_27 (or .def_26 .def_23))) (let ((.def_28 (and .def_27 .def_21))) (let ((.def_29 (not .def_28))) (let ((.def_30 (and .def_29 .def_14))) (let ((.def_31 (and A24 A19))) (let ((.def_32 (not .def_31))) (let ((.def_33 (and A10 A15))) (let ((.def_34 (or .def_33 .def_32))) (let ((.def_35 (not .def_34))) (let ((.def_36 (not A1))) (let ((.def_37 (and .def_36 A24))) (let ((.def_38 (not .def_37))) (let ((.def_39 (not A10))) (let ((.def_40 (and .def_39 .def_0))) (let ((.def_41 (not .def_40))) (let ((.def_42 (or .def_41 .def_38))) (let ((.def_43 (not .def_42))) (let ((.def_44 (and .def_43 .def_35))) (let ((.def_45 (not .def_44))) (let ((.def_46 (and .def_3 .def_22))) (let ((.def_47 (not .def_46))) (let ((.def_48 (not A11))) (let ((.def_49 (or .def_48 A15))) (let ((.def_50 (not .def_49))) (let ((.def_51 (and .def_50 .def_47))) (let ((.def_52 (not .def_51))) (let ((.def_53 (and .def_19 A10))) (let ((.def_54 (not .def_53))) (let ((.def_55 (not A9))) (let ((.def_56 (and A10 .def_55))) (let ((.def_57 (not .def_56))) (let ((.def_58 (and .def_57 .def_54))) (let ((.def_59 (not .def_58))) (let ((.def_60 (= .def_59 .def_52))) (let ((.def_61 (not .def_60))) (let ((.def_62 (or .def_61 .def_45))) (let ((.def_63 (not .def_62))) (let ((.def_64 (and .def_63 .def_30))) (let ((.def_65 (not .def_64))) (let ((.def_66 (not A14))) (let ((.def_67 (= A0 .def_66))) (let ((.def_68 (not .def_67))) (let ((.def_69 (or .def_15 .def_24))) (let ((.def_70 (not .def_69))) (let ((.def_71 (and .def_70 .def_68))) (let ((.def_72 (not .def_71))) (let ((.def_73 (= A13 A18))) (let ((.def_74 (not .def_73))) (let ((.def_75 (not A13))) (let ((.def_76 (not A21))) (let ((.def_77 (or .def_76 .def_75))) (let ((.def_78 (not .def_77))) (let ((.def_79 (and .def_78 .def_74))) (let ((.def_80 (not .def_79))) (let ((.def_81 (or .def_80 .def_72))) (let ((.def_82 (not .def_81))) (let ((.def_83 (= .def_0 A18))) (let ((.def_84 (not .def_83))) (let ((.def_85 (or A23 A8))) (let ((.def_86 (and .def_85 .def_84))) (let ((.def_87 (not .def_86))) (let ((.def_88 (= A1 A24))) (let ((.def_89 (or A11 A16))) (let ((.def_90 (not .def_89))) (let ((.def_91 (and .def_90 .def_88))) (let ((.def_92 (or .def_91 .def_87))) (let ((.def_93 (or .def_92 .def_82))) (let ((.def_94 (or A20 .def_18))) (let ((.def_95 (and .def_55 A3))) (let ((.def_96 (not .def_95))) (let ((.def_97 (= .def_96 .def_94))) (let ((.def_98 (not .def_97))) (let ((.def_99 (not A8))) (let ((.def_100 (and .def_99 A17))) (let ((.def_101 (not .def_100))) (let ((.def_102 (and .def_3 A3))) (let ((.def_103 (or .def_102 .def_101))) (let ((.def_104 (and .def_103 .def_98))) (let ((.def_105 (not A5))) (let ((.def_106 (and A22 .def_105))) (let ((.def_107 (not .def_106))) (let ((.def_108 (= .def_99 A9))) (let ((.def_109 (not .def_108))) (let ((.def_110 (and .def_109 .def_107))) (let ((.def_111 (not A24))) (let ((.def_112 (and .def_111 .def_75))) (let ((.def_113 (not .def_112))) (let ((.def_114 (or .def_76 A1))) (let ((.def_115 (or .def_114 .def_113))) (let ((.def_116 (or .def_115 .def_110))) (let ((.def_117 (and .def_116 .def_104))) (let ((.def_118 (not .def_117))) (let ((.def_119 (and .def_118 .def_93))) (let ((.def_120 (= .def_119 .def_65))) .def_120))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
