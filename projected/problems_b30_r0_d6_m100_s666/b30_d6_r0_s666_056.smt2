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
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
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
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A21))) (let ((.def_1 (not A28))) (let ((.def_2 (and .def_1 .def_0))) (let ((.def_3 (not .def_2))) (let ((.def_4 (and A6 A28))) (let ((.def_5 (or .def_4 .def_3))) (let ((.def_6 (not .def_5))) (let ((.def_7 (not A8))) (let ((.def_8 (and .def_7 A29))) (let ((.def_9 (not A23))) (let ((.def_10 (not A20))) (let ((.def_11 (or .def_10 .def_9))) (let ((.def_12 (and .def_11 .def_8))) (let ((.def_13 (not .def_12))) (let ((.def_14 (or .def_13 .def_6))) (let ((.def_15 (not .def_14))) (let ((.def_16 (not A16))) (let ((.def_17 (and A1 .def_16))) (let ((.def_18 (not A1))) (let ((.def_19 (or .def_18 A18))) (let ((.def_20 (and .def_19 .def_17))) (let ((.def_21 (not .def_20))) (let ((.def_22 (or A20 .def_9))) (let ((.def_23 (not .def_22))) (let ((.def_24 (or A13 A13))) (let ((.def_25 (not .def_24))) (let ((.def_26 (or .def_25 .def_23))) (let ((.def_27 (or .def_26 .def_21))) (let ((.def_28 (not .def_27))) (let ((.def_29 (and .def_28 .def_15))) (let ((.def_30 (not .def_29))) (let ((.def_31 (and .def_0 .def_10))) (let ((.def_32 (not .def_31))) (let ((.def_33 (= A7 A18))) (let ((.def_34 (and .def_33 .def_32))) (let ((.def_35 (not .def_34))) (let ((.def_36 (and A12 A0))) (let ((.def_37 (not .def_36))) (let ((.def_38 (not A0))) (let ((.def_39 (or A24 .def_38))) (let ((.def_40 (not .def_39))) (let ((.def_41 (and .def_40 .def_37))) (let ((.def_42 (not .def_41))) (let ((.def_43 (or .def_42 .def_35))) (let ((.def_44 (not .def_43))) (let ((.def_45 (or .def_38 A21))) (let ((.def_46 (not .def_45))) (let ((.def_47 (not A7))) (let ((.def_48 (not A5))) (let ((.def_49 (and .def_48 .def_47))) (let ((.def_50 (not .def_49))) (let ((.def_51 (and .def_50 .def_46))) (let ((.def_52 (and A2 A24))) (let ((.def_53 (not .def_52))) (let ((.def_54 (and A28 .def_10))) (let ((.def_55 (not .def_54))) (let ((.def_56 (or .def_55 .def_53))) (let ((.def_57 (not .def_56))) (let ((.def_58 (or .def_57 .def_51))) (let ((.def_59 (not .def_58))) (let ((.def_60 (or .def_59 .def_44))) (let ((.def_61 (and .def_60 .def_30))) (let ((.def_62 (not A27))) (let ((.def_63 (and .def_62 .def_48))) (let ((.def_64 (not A12))) (let ((.def_65 (and A24 .def_64))) (let ((.def_66 (= .def_65 .def_63))) (let ((.def_67 (not .def_66))) (let ((.def_68 (not A2))) (let ((.def_69 (or A16 .def_68))) (let ((.def_70 (not .def_69))) (let ((.def_71 (or A29 A29))) (let ((.def_72 (or .def_71 .def_70))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and .def_73 .def_67))) (let ((.def_75 (not A17))) (let ((.def_76 (not A3))) (let ((.def_77 (= .def_76 .def_75))) (let ((.def_78 (and A25 A17))) (let ((.def_79 (not .def_78))) (let ((.def_80 (or .def_79 .def_77))) (let ((.def_81 (or A11 A27))) (let ((.def_82 (not .def_81))) (let ((.def_83 (not A11))) (let ((.def_84 (and .def_38 .def_83))) (let ((.def_85 (or .def_84 .def_82))) (let ((.def_86 (not .def_85))) (let ((.def_87 (or .def_86 .def_80))) (let ((.def_88 (not .def_87))) (let ((.def_89 (= .def_88 .def_74))) (let ((.def_90 (not .def_89))) (let ((.def_91 (not A24))) (let ((.def_92 (or .def_76 .def_91))) (let ((.def_93 (not .def_92))) (let ((.def_94 (not A10))) (let ((.def_95 (and .def_94 A22))) (let ((.def_96 (not .def_95))) (let ((.def_97 (or .def_96 .def_93))) (let ((.def_98 (or A19 .def_94))) (let ((.def_99 (not .def_98))) (let ((.def_100 (or A5 .def_0))) (let ((.def_101 (or .def_100 .def_99))) (let ((.def_102 (not .def_101))) (let ((.def_103 (or .def_102 .def_97))) (let ((.def_104 (not A19))) (let ((.def_105 (and A20 .def_104))) (let ((.def_106 (and A0 A29))) (let ((.def_107 (and .def_106 .def_105))) (let ((.def_108 (not .def_107))) (let ((.def_109 (and A4 A25))) (let ((.def_110 (not .def_109))) (let ((.def_111 (or .def_0 A17))) (let ((.def_112 (or .def_111 .def_110))) (let ((.def_113 (or .def_112 .def_108))) (let ((.def_114 (not .def_113))) (let ((.def_115 (or .def_114 .def_103))) (let ((.def_116 (not .def_115))) (let ((.def_117 (and .def_116 .def_90))) (let ((.def_118 (and .def_117 .def_61))) .def_118))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
