(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
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
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (not A2))) (let ((.def_1 (not A0))) (let ((.def_2 (and .def_1 .def_0))) (let ((.def_3 (not A13))) (let ((.def_4 (and A13 .def_3))) (let ((.def_5 (not .def_4))) (let ((.def_6 (and .def_5 .def_2))) (let ((.def_7 (not A1))) (let ((.def_8 (and A19 .def_7))) (let ((.def_9 (not A4))) (let ((.def_10 (not A23))) (let ((.def_11 (and .def_10 .def_9))) (let ((.def_12 (not .def_11))) (let ((.def_13 (or .def_12 .def_8))) (let ((.def_14 (not .def_13))) (let ((.def_15 (or .def_14 .def_6))) (let ((.def_16 (not .def_15))) (let ((.def_17 (or A21 A12))) (let ((.def_18 (not A8))) (let ((.def_19 (or .def_18 A9))) (let ((.def_20 (or .def_19 .def_17))) (let ((.def_21 (not .def_20))) (let ((.def_22 (not A20))) (let ((.def_23 (and A19 .def_22))) (let ((.def_24 (or A14 A14))) (let ((.def_25 (not .def_24))) (let ((.def_26 (or .def_25 .def_23))) (let ((.def_27 (not .def_26))) (let ((.def_28 (and .def_27 .def_21))) (let ((.def_29 (not .def_28))) (let ((.def_30 (or .def_29 .def_16))) (let ((.def_31 (and A10 .def_10))) (let ((.def_32 (not .def_31))) (let ((.def_33 (= .def_10 A21))) (let ((.def_34 (or .def_33 .def_32))) (let ((.def_35 (not A7))) (let ((.def_36 (and .def_35 .def_7))) (let ((.def_37 (not A14))) (let ((.def_38 (and .def_37 A4))) (let ((.def_39 (and .def_38 .def_36))) (let ((.def_40 (not .def_39))) (let ((.def_41 (or .def_40 .def_34))) (let ((.def_42 (not A3))) (let ((.def_43 (or A14 .def_42))) (let ((.def_44 (or A19 A12))) (let ((.def_45 (not .def_44))) (let ((.def_46 (and .def_45 .def_43))) (let ((.def_47 (not A22))) (let ((.def_48 (or .def_47 A24))) (let ((.def_49 (and .def_3 A6))) (let ((.def_50 (not .def_49))) (let ((.def_51 (and .def_50 .def_48))) (let ((.def_52 (or .def_51 .def_46))) (let ((.def_53 (not .def_52))) (let ((.def_54 (or .def_53 .def_41))) (let ((.def_55 (not .def_54))) (let ((.def_56 (and .def_55 .def_30))) (let ((.def_57 (not .def_56))) (let ((.def_58 (and .def_3 A21))) (let ((.def_59 (or A9 A22))) (let ((.def_60 (not .def_59))) (let ((.def_61 (or .def_60 .def_58))) (let ((.def_62 (or A2 A17))) (let ((.def_63 (not A11))) (let ((.def_64 (or .def_63 A1))) (let ((.def_65 (and .def_64 .def_62))) (let ((.def_66 (and .def_65 .def_61))) (let ((.def_67 (not .def_66))) (let ((.def_68 (= A12 A18))) (let ((.def_69 (not .def_68))) (let ((.def_70 (and A24 A9))) (let ((.def_71 (not .def_70))) (let ((.def_72 (and .def_71 .def_69))) (let ((.def_73 (not .def_72))) (let ((.def_74 (or A18 A14))) (let ((.def_75 (not .def_74))) (let ((.def_76 (not A5))) (let ((.def_77 (not A16))) (let ((.def_78 (and .def_77 .def_76))) (let ((.def_79 (and .def_78 .def_75))) (let ((.def_80 (not .def_79))) (let ((.def_81 (or .def_80 .def_73))) (let ((.def_82 (not .def_81))) (let ((.def_83 (= .def_82 .def_67))) (let ((.def_84 (not .def_83))) (let ((.def_85 (and .def_0 A14))) (let ((.def_86 (not A9))) (let ((.def_87 (or .def_86 A2))) (let ((.def_88 (or .def_87 .def_85))) (let ((.def_89 (or .def_7 A24))) (let ((.def_90 (not .def_89))) (let ((.def_91 (or A18 A5))) (let ((.def_92 (not .def_91))) (let ((.def_93 (or .def_92 .def_90))) (let ((.def_94 (and .def_93 .def_88))) (let ((.def_95 (or .def_3 .def_1))) (let ((.def_96 (= .def_18 A8))) (let ((.def_97 (not .def_96))) (let ((.def_98 (and .def_97 .def_95))) (let ((.def_99 (not .def_98))) (let ((.def_100 (or A10 A3))) (let ((.def_101 (not .def_100))) (let ((.def_102 (and .def_10 A22))) (let ((.def_103 (not .def_102))) (let ((.def_104 (or .def_103 .def_101))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and .def_105 .def_99))) (let ((.def_107 (not .def_106))) (let ((.def_108 (= .def_107 .def_94))) (let ((.def_109 (not .def_108))) (let ((.def_110 (and .def_109 .def_84))) (let ((.def_111 (not .def_110))) (let ((.def_112 (and .def_111 .def_57))) (let ((.def_113 (not .def_112))) .def_113)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
