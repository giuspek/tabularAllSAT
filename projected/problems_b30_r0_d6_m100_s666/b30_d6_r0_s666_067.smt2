(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
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
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A19))) (let ((.def_1 (or .def_0 A2))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A26))) (let ((.def_4 (= .def_3 A0))) (let ((.def_5 (or .def_4 .def_2))) (let ((.def_6 (not .def_5))) (let ((.def_7 (not A17))) (let ((.def_8 (not A12))) (let ((.def_9 (and .def_8 .def_7))) (let ((.def_10 (not A11))) (let ((.def_11 (and A29 .def_10))) (let ((.def_12 (not .def_11))) (let ((.def_13 (or .def_12 .def_9))) (let ((.def_14 (or .def_13 .def_6))) (let ((.def_15 (not A7))) (let ((.def_16 (not A21))) (let ((.def_17 (or .def_16 .def_15))) (let ((.def_18 (not A1))) (let ((.def_19 (not A5))) (let ((.def_20 (or .def_19 .def_18))) (let ((.def_21 (not .def_20))) (let ((.def_22 (or .def_21 .def_17))) (let ((.def_23 (not A18))) (let ((.def_24 (and .def_23 A19))) (let ((.def_25 (not .def_24))) (let ((.def_26 (and A2 A17))) (let ((.def_27 (and .def_26 .def_25))) (let ((.def_28 (or .def_27 .def_22))) (let ((.def_29 (not .def_28))) (let ((.def_30 (= .def_29 .def_14))) (let ((.def_31 (not .def_30))) (let ((.def_32 (not A15))) (let ((.def_33 (not A23))) (let ((.def_34 (or .def_33 .def_32))) (let ((.def_35 (not .def_34))) (let ((.def_36 (not A27))) (let ((.def_37 (not A10))) (let ((.def_38 (or .def_37 .def_36))) (let ((.def_39 (or .def_38 .def_35))) (let ((.def_40 (or A0 A0))) (let ((.def_41 (or .def_16 A24))) (let ((.def_42 (= .def_41 .def_40))) (let ((.def_43 (not .def_42))) (let ((.def_44 (or .def_43 .def_39))) (let ((.def_45 (not A6))) (let ((.def_46 (or .def_45 A12))) (let ((.def_47 (or .def_33 A13))) (let ((.def_48 (or .def_47 .def_46))) (let ((.def_49 (not .def_48))) (let ((.def_50 (not A9))) (let ((.def_51 (or A18 .def_50))) (let ((.def_52 (not .def_51))) (let ((.def_53 (and .def_50 A27))) (let ((.def_54 (= .def_53 .def_52))) (let ((.def_55 (or .def_54 .def_49))) (let ((.def_56 (= .def_55 .def_44))) (let ((.def_57 (not .def_56))) (let ((.def_58 (or .def_57 .def_31))) (let ((.def_59 (not .def_58))) (let ((.def_60 (not A16))) (let ((.def_61 (or A17 .def_60))) (let ((.def_62 (not A28))) (let ((.def_63 (or .def_62 A28))) (let ((.def_64 (not .def_63))) (let ((.def_65 (and .def_64 .def_61))) (let ((.def_66 (and A5 A8))) (let ((.def_67 (not .def_66))) (let ((.def_68 (not A0))) (let ((.def_69 (and A22 .def_68))) (let ((.def_70 (not .def_69))) (let ((.def_71 (or .def_70 .def_67))) (let ((.def_72 (and .def_71 .def_65))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and A15 A8))) (let ((.def_75 (not A13))) (let ((.def_76 (and A21 .def_75))) (let ((.def_77 (and .def_76 .def_74))) (let ((.def_78 (not .def_77))) (let ((.def_79 (or A24 .def_33))) (let ((.def_80 (or .def_16 A21))) (let ((.def_81 (not .def_80))) (let ((.def_82 (= .def_81 .def_79))) (let ((.def_83 (not .def_82))) (let ((.def_84 (= .def_83 .def_78))) (let ((.def_85 (not .def_84))) (let ((.def_86 (or .def_85 .def_73))) (let ((.def_87 (not .def_86))) (let ((.def_88 (or .def_68 .def_15))) (let ((.def_89 (and A24 .def_32))) (let ((.def_90 (and .def_89 .def_88))) (let ((.def_91 (not .def_90))) (let ((.def_92 (and A28 .def_19))) (let ((.def_93 (not .def_92))) (let ((.def_94 (and .def_10 .def_23))) (let ((.def_95 (and .def_94 .def_93))) (let ((.def_96 (not .def_95))) (let ((.def_97 (or .def_96 .def_91))) (let ((.def_98 (not A2))) (let ((.def_99 (not A24))) (let ((.def_100 (and .def_99 .def_98))) (let ((.def_101 (not .def_100))) (let ((.def_102 (or .def_45 A8))) (let ((.def_103 (or .def_102 .def_101))) (let ((.def_104 (= A23 A26))) (let ((.def_105 (and .def_7 .def_3))) (let ((.def_106 (or .def_105 .def_104))) (let ((.def_107 (not .def_106))) (let ((.def_108 (or .def_107 .def_103))) (let ((.def_109 (or .def_108 .def_97))) (let ((.def_110 (not .def_109))) (let ((.def_111 (and .def_110 .def_87))) (let ((.def_112 (and .def_111 .def_59))) .def_112))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
