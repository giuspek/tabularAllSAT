(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
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
(assert (let ((.def_0 (not A9))) (let ((.def_1 (= .def_0 A5))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A26))) (let ((.def_4 (and A25 .def_3))) (let ((.def_5 (not .def_4))) (let ((.def_6 (or .def_5 .def_2))) (let ((.def_7 (not A2))) (let ((.def_8 (not A12))) (let ((.def_9 (and .def_8 .def_7))) (let ((.def_10 (not .def_9))) (let ((.def_11 (or A29 A19))) (let ((.def_12 (or .def_11 .def_10))) (let ((.def_13 (= .def_12 .def_6))) (let ((.def_14 (not A24))) (let ((.def_15 (not A22))) (let ((.def_16 (or .def_15 .def_14))) (let ((.def_17 (and A2 A26))) (let ((.def_18 (not .def_17))) (let ((.def_19 (or .def_18 .def_16))) (let ((.def_20 (not .def_19))) (let ((.def_21 (not A18))) (let ((.def_22 (or A17 .def_21))) (let ((.def_23 (and A28 A18))) (let ((.def_24 (not .def_23))) (let ((.def_25 (and .def_24 .def_22))) (let ((.def_26 (not .def_25))) (let ((.def_27 (and .def_26 .def_20))) (let ((.def_28 (or .def_27 .def_13))) (let ((.def_29 (not .def_28))) (let ((.def_30 (and A13 .def_21))) (let ((.def_31 (not .def_30))) (let ((.def_32 (not A28))) (let ((.def_33 (not A25))) (let ((.def_34 (or .def_33 .def_32))) (let ((.def_35 (not .def_34))) (let ((.def_36 (= .def_35 .def_31))) (let ((.def_37 (not .def_36))) (let ((.def_38 (not A27))) (let ((.def_39 (or .def_38 A19))) (let ((.def_40 (not .def_39))) (let ((.def_41 (not A16))) (let ((.def_42 (not A4))) (let ((.def_43 (or .def_42 .def_41))) (let ((.def_44 (not .def_43))) (let ((.def_45 (and .def_44 .def_40))) (let ((.def_46 (or .def_45 .def_37))) (let ((.def_47 (and A16 A3))) (let ((.def_48 (not A11))) (let ((.def_49 (and .def_48 A0))) (let ((.def_50 (not .def_49))) (let ((.def_51 (= .def_50 .def_47))) (let ((.def_52 (not A1))) (let ((.def_53 (or .def_52 .def_38))) (let ((.def_54 (not .def_53))) (let ((.def_55 (or .def_41 .def_3))) (let ((.def_56 (or .def_55 .def_54))) (let ((.def_57 (or .def_56 .def_51))) (let ((.def_58 (not .def_57))) (let ((.def_59 (= .def_58 .def_46))) (let ((.def_60 (or .def_59 .def_29))) (let ((.def_61 (not A29))) (let ((.def_62 (= A16 .def_61))) (let ((.def_63 (not A19))) (let ((.def_64 (and .def_63 A29))) (let ((.def_65 (or .def_64 .def_62))) (let ((.def_66 (not .def_65))) (let ((.def_67 (or A25 A25))) (let ((.def_68 (and A26 A8))) (let ((.def_69 (and .def_68 .def_67))) (let ((.def_70 (not .def_69))) (let ((.def_71 (and .def_70 .def_66))) (let ((.def_72 (or .def_48 A22))) (let ((.def_73 (not A0))) (let ((.def_74 (not A17))) (let ((.def_75 (or .def_74 .def_73))) (let ((.def_76 (not .def_75))) (let ((.def_77 (and .def_76 .def_72))) (let ((.def_78 (or A10 A23))) (let ((.def_79 (not A20))) (let ((.def_80 (= .def_15 .def_79))) (let ((.def_81 (not .def_80))) (let ((.def_82 (= .def_81 .def_78))) (let ((.def_83 (not .def_82))) (let ((.def_84 (and .def_83 .def_77))) (let ((.def_85 (or .def_84 .def_71))) (let ((.def_86 (not .def_85))) (let ((.def_87 (and A12 A15))) (let ((.def_88 (and .def_32 A16))) (let ((.def_89 (not .def_88))) (let ((.def_90 (or .def_89 .def_87))) (let ((.def_91 (or A9 .def_61))) (let ((.def_92 (not .def_91))) (let ((.def_93 (= A22 A11))) (let ((.def_94 (or .def_93 .def_92))) (let ((.def_95 (or .def_94 .def_90))) (let ((.def_96 (not A3))) (let ((.def_97 (or .def_96 A20))) (let ((.def_98 (not A8))) (let ((.def_99 (= .def_98 A10))) (let ((.def_100 (not .def_99))) (let ((.def_101 (and .def_100 .def_97))) (let ((.def_102 (= .def_3 A26))) (let ((.def_103 (not .def_102))) (let ((.def_104 (and A1 A26))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and .def_105 .def_103))) (let ((.def_107 (not .def_106))) (let ((.def_108 (or .def_107 .def_101))) (let ((.def_109 (not .def_108))) (let ((.def_110 (and .def_109 .def_95))) (let ((.def_111 (not .def_110))) (let ((.def_112 (or .def_111 .def_86))) (let ((.def_113 (not .def_112))) (let ((.def_114 (or .def_113 .def_60))) (let ((.def_115 (not .def_114))) .def_115)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
