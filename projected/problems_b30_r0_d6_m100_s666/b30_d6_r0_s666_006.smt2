(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
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
(declare-fun A24 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A0))) (let ((.def_1 (and .def_0 A21))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A10))) (let ((.def_4 (= A8 .def_3))) (let ((.def_5 (or .def_4 .def_2))) (let ((.def_6 (not .def_5))) (let ((.def_7 (and A17 A1))) (let ((.def_8 (not .def_7))) (let ((.def_9 (not A8))) (let ((.def_10 (and .def_0 .def_9))) (let ((.def_11 (not .def_10))) (let ((.def_12 (or .def_11 .def_8))) (let ((.def_13 (or .def_12 .def_6))) (let ((.def_14 (= A29 .def_9))) (let ((.def_15 (not A2))) (let ((.def_16 (not A6))) (let ((.def_17 (or .def_16 .def_15))) (let ((.def_18 (and .def_17 .def_14))) (let ((.def_19 (not .def_18))) (let ((.def_20 (not A27))) (let ((.def_21 (or .def_3 .def_20))) (let ((.def_22 (not A12))) (let ((.def_23 (or A9 .def_22))) (let ((.def_24 (or .def_23 .def_21))) (let ((.def_25 (not .def_24))) (let ((.def_26 (or .def_25 .def_19))) (let ((.def_27 (and .def_26 .def_13))) (let ((.def_28 (not .def_27))) (let ((.def_29 (and A2 .def_20))) (let ((.def_30 (and A26 A14))) (let ((.def_31 (not .def_30))) (let ((.def_32 (and .def_31 .def_29))) (let ((.def_33 (not .def_32))) (let ((.def_34 (or A12 .def_15))) (let ((.def_35 (not .def_34))) (let ((.def_36 (and A22 A13))) (let ((.def_37 (or .def_36 .def_35))) (let ((.def_38 (not .def_37))) (let ((.def_39 (and .def_38 .def_33))) (let ((.def_40 (not A9))) (let ((.def_41 (not A21))) (let ((.def_42 (and .def_41 .def_40))) (let ((.def_43 (not .def_42))) (let ((.def_44 (or .def_20 A18))) (let ((.def_45 (or .def_44 .def_43))) (let ((.def_46 (not .def_45))) (let ((.def_47 (or A20 A19))) (let ((.def_48 (not .def_47))) (let ((.def_49 (or .def_0 A8))) (let ((.def_50 (or .def_49 .def_48))) (let ((.def_51 (or .def_50 .def_46))) (let ((.def_52 (not .def_51))) (let ((.def_53 (and .def_52 .def_39))) (let ((.def_54 (or .def_53 .def_28))) (let ((.def_55 (not A28))) (let ((.def_56 (not A16))) (let ((.def_57 (or .def_56 .def_55))) (let ((.def_58 (not .def_57))) (let ((.def_59 (= A6 A7))) (let ((.def_60 (and .def_59 .def_58))) (let ((.def_61 (and A29 A16))) (let ((.def_62 (not A13))) (let ((.def_63 (= .def_62 A28))) (let ((.def_64 (or .def_63 .def_61))) (let ((.def_65 (or .def_64 .def_60))) (let ((.def_66 (or .def_22 .def_16))) (let ((.def_67 (not .def_66))) (let ((.def_68 (or .def_56 A28))) (let ((.def_69 (= .def_68 .def_67))) (let ((.def_70 (not A22))) (let ((.def_71 (or .def_70 A4))) (let ((.def_72 (and A4 A18))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and .def_73 .def_71))) (let ((.def_75 (or .def_74 .def_69))) (let ((.def_76 (not .def_75))) (let ((.def_77 (and .def_76 .def_65))) (let ((.def_78 (not A24))) (let ((.def_79 (and A16 .def_78))) (let ((.def_80 (not .def_79))) (let ((.def_81 (or .def_62 .def_41))) (let ((.def_82 (not .def_81))) (let ((.def_83 (and .def_82 .def_80))) (let ((.def_84 (not .def_83))) (let ((.def_85 (or .def_78 .def_20))) (let ((.def_86 (not A11))) (let ((.def_87 (and .def_22 .def_86))) (let ((.def_88 (or .def_87 .def_85))) (let ((.def_89 (and .def_88 .def_84))) (let ((.def_90 (not A29))) (let ((.def_91 (or .def_40 .def_90))) (let ((.def_92 (not A26))) (let ((.def_93 (or .def_92 .def_3))) (let ((.def_94 (or .def_93 .def_91))) (let ((.def_95 (not .def_94))) (let ((.def_96 (not A19))) (let ((.def_97 (or .def_96 .def_0))) (let ((.def_98 (not .def_97))) (let ((.def_99 (not A5))) (let ((.def_100 (= .def_99 A28))) (let ((.def_101 (or .def_100 .def_98))) (let ((.def_102 (not .def_101))) (let ((.def_103 (and .def_102 .def_95))) (let ((.def_104 (and .def_103 .def_89))) (let ((.def_105 (or .def_104 .def_77))) (let ((.def_106 (and .def_105 .def_54))) (let ((.def_107 (not .def_106))) .def_107)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
