(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
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
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (not A24))) (let ((.def_1 (not A14))) (let ((.def_2 (= .def_1 .def_0))) (let ((.def_3 (not .def_2))) (let ((.def_4 (not A1))) (let ((.def_5 (and A19 .def_4))) (let ((.def_6 (and .def_5 .def_3))) (let ((.def_7 (not A15))) (let ((.def_8 (not A8))) (let ((.def_9 (or .def_8 .def_7))) (let ((.def_10 (not .def_9))) (let ((.def_11 (not A23))) (let ((.def_12 (or .def_11 A2))) (let ((.def_13 (not .def_12))) (let ((.def_14 (and .def_13 .def_10))) (let ((.def_15 (not .def_14))) (let ((.def_16 (= .def_15 .def_6))) (let ((.def_17 (not .def_16))) (let ((.def_18 (not A7))) (let ((.def_19 (or .def_18 .def_4))) (let ((.def_20 (not .def_19))) (let ((.def_21 (or A8 .def_1))) (let ((.def_22 (and .def_21 .def_20))) (let ((.def_23 (not .def_22))) (let ((.def_24 (and A2 A7))) (let ((.def_25 (not .def_24))) (let ((.def_26 (or A12 A9))) (let ((.def_27 (not .def_26))) (let ((.def_28 (or .def_27 .def_25))) (let ((.def_29 (and .def_28 .def_23))) (let ((.def_30 (= .def_29 .def_17))) (let ((.def_31 (not .def_30))) (let ((.def_32 (not A9))) (let ((.def_33 (not A13))) (let ((.def_34 (or .def_33 .def_32))) (let ((.def_35 (= .def_4 A22))) (let ((.def_36 (not .def_35))) (let ((.def_37 (= .def_36 .def_34))) (let ((.def_38 (not .def_37))) (let ((.def_39 (and A14 .def_18))) (let ((.def_40 (or A16 .def_1))) (let ((.def_41 (= .def_40 .def_39))) (let ((.def_42 (= .def_41 .def_38))) (let ((.def_43 (or A0 A2))) (let ((.def_44 (not .def_43))) (let ((.def_45 (not A2))) (let ((.def_46 (and .def_45 A24))) (let ((.def_47 (not .def_46))) (let ((.def_48 (and .def_47 .def_44))) (let ((.def_49 (not .def_48))) (let ((.def_50 (not A3))) (let ((.def_51 (not A21))) (let ((.def_52 (and .def_51 .def_50))) (let ((.def_53 (or A22 .def_45))) (let ((.def_54 (not .def_53))) (let ((.def_55 (or .def_54 .def_52))) (let ((.def_56 (not .def_55))) (let ((.def_57 (and .def_56 .def_49))) (let ((.def_58 (not .def_57))) (let ((.def_59 (or .def_58 .def_42))) (let ((.def_60 (and .def_59 .def_31))) (let ((.def_61 (not .def_60))) (let ((.def_62 (or A0 A24))) (let ((.def_63 (not .def_62))) (let ((.def_64 (or A20 A2))) (let ((.def_65 (and .def_64 .def_63))) (let ((.def_66 (or A2 A4))) (let ((.def_67 (= .def_45 A3))) (let ((.def_68 (= .def_67 .def_66))) (let ((.def_69 (not .def_68))) (let ((.def_70 (or .def_69 .def_65))) (let ((.def_71 (not .def_70))) (let ((.def_72 (and A4 A1))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and A19 .def_1))) (let ((.def_75 (not .def_74))) (let ((.def_76 (and .def_75 .def_73))) (let ((.def_77 (not .def_76))) (let ((.def_78 (not A19))) (let ((.def_79 (or .def_78 A4))) (let ((.def_80 (not .def_79))) (let ((.def_81 (not A17))) (let ((.def_82 (or .def_81 A0))) (let ((.def_83 (and .def_82 .def_80))) (let ((.def_84 (not .def_83))) (let ((.def_85 (and .def_84 .def_77))) (let ((.def_86 (and .def_85 .def_71))) (let ((.def_87 (not .def_86))) (let ((.def_88 (or A4 A7))) (let ((.def_89 (or .def_51 A1))) (let ((.def_90 (and .def_89 .def_88))) (let ((.def_91 (not .def_90))) (let ((.def_92 (and A8 A10))) (let ((.def_93 (not .def_92))) (let ((.def_94 (and .def_18 A10))) (let ((.def_95 (not .def_94))) (let ((.def_96 (and .def_95 .def_93))) (let ((.def_97 (and .def_96 .def_91))) (let ((.def_98 (and A13 .def_7))) (let ((.def_99 (not A20))) (let ((.def_100 (and .def_99 .def_50))) (let ((.def_101 (not .def_100))) (let ((.def_102 (or .def_101 .def_98))) (let ((.def_103 (not A11))) (let ((.def_104 (= A2 .def_103))) (let ((.def_105 (not .def_104))) (let ((.def_106 (not A22))) (let ((.def_107 (and .def_99 .def_106))) (let ((.def_108 (and .def_107 .def_105))) (let ((.def_109 (or .def_108 .def_102))) (let ((.def_110 (= .def_109 .def_97))) (let ((.def_111 (or .def_110 .def_87))) (let ((.def_112 (not .def_111))) (let ((.def_113 (or .def_112 .def_61))) (let ((.def_114 (not .def_113))) .def_114))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
