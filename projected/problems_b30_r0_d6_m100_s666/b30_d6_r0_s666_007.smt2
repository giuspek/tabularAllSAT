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
(declare-fun A21 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(declare-fun A25 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A19))) (let ((.def_1 (not A0))) (let ((.def_2 (or .def_1 .def_0))) (let ((.def_3 (not A24))) (let ((.def_4 (or A9 .def_3))) (let ((.def_5 (not .def_4))) (let ((.def_6 (or .def_5 .def_2))) (let ((.def_7 (not .def_6))) (let ((.def_8 (= A14 A26))) (let ((.def_9 (or A10 A18))) (let ((.def_10 (and .def_9 .def_8))) (let ((.def_11 (or .def_10 .def_7))) (let ((.def_12 (not .def_11))) (let ((.def_13 (and .def_1 A0))) (let ((.def_14 (not A16))) (let ((.def_15 (= .def_14 A17))) (let ((.def_16 (and .def_15 .def_13))) (let ((.def_17 (not .def_16))) (let ((.def_18 (or A27 A16))) (let ((.def_19 (not .def_18))) (let ((.def_20 (not A25))) (let ((.def_21 (and A14 .def_20))) (let ((.def_22 (or .def_21 .def_19))) (let ((.def_23 (or .def_22 .def_17))) (let ((.def_24 (and .def_23 .def_12))) (let ((.def_25 (not A28))) (let ((.def_26 (not A11))) (let ((.def_27 (or .def_26 .def_25))) (let ((.def_28 (not A1))) (let ((.def_29 (or A25 .def_28))) (let ((.def_30 (or .def_29 .def_27))) (let ((.def_31 (or A1 A10))) (let ((.def_32 (not .def_31))) (let ((.def_33 (and A0 A26))) (let ((.def_34 (not .def_33))) (let ((.def_35 (and .def_34 .def_32))) (let ((.def_36 (and .def_35 .def_30))) (let ((.def_37 (not .def_36))) (let ((.def_38 (not A13))) (let ((.def_39 (or .def_28 .def_38))) (let ((.def_40 (not A2))) (let ((.def_41 (and .def_40 A23))) (let ((.def_42 (or .def_41 .def_39))) (let ((.def_43 (not .def_42))) (let ((.def_44 (not A4))) (let ((.def_45 (= .def_44 A7))) (let ((.def_46 (and A24 A19))) (let ((.def_47 (not .def_46))) (let ((.def_48 (or .def_47 .def_45))) (let ((.def_49 (not .def_48))) (let ((.def_50 (or .def_49 .def_43))) (let ((.def_51 (not .def_50))) (let ((.def_52 (= .def_51 .def_37))) (let ((.def_53 (not .def_52))) (let ((.def_54 (and .def_53 .def_24))) (let ((.def_55 (or A7 A5))) (let ((.def_56 (not .def_55))) (let ((.def_57 (not A8))) (let ((.def_58 (not A29))) (let ((.def_59 (and .def_58 .def_57))) (let ((.def_60 (or .def_59 .def_56))) (let ((.def_61 (not .def_60))) (let ((.def_62 (= .def_14 A26))) (let ((.def_63 (not .def_62))) (let ((.def_64 (not A12))) (let ((.def_65 (or A2 .def_64))) (let ((.def_66 (= .def_65 .def_63))) (let ((.def_67 (and .def_66 .def_61))) (let ((.def_68 (= .def_25 A29))) (let ((.def_69 (not .def_68))) (let ((.def_70 (not A15))) (let ((.def_71 (and .def_70 A28))) (let ((.def_72 (or .def_71 .def_69))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and .def_44 A2))) (let ((.def_75 (not .def_74))) (let ((.def_76 (not A21))) (let ((.def_77 (not A5))) (let ((.def_78 (and .def_77 .def_76))) (let ((.def_79 (or .def_78 .def_75))) (let ((.def_80 (not .def_79))) (let ((.def_81 (or .def_80 .def_73))) (let ((.def_82 (and .def_81 .def_67))) (let ((.def_83 (not .def_82))) (let ((.def_84 (not A17))) (let ((.def_85 (or A24 .def_84))) (let ((.def_86 (not .def_85))) (let ((.def_87 (and .def_40 .def_14))) (let ((.def_88 (not .def_87))) (let ((.def_89 (and .def_88 .def_86))) (let ((.def_90 (not .def_89))) (let ((.def_91 (or A8 A2))) (let ((.def_92 (not .def_91))) (let ((.def_93 (not A14))) (let ((.def_94 (not A3))) (let ((.def_95 (or .def_94 .def_93))) (let ((.def_96 (or .def_95 .def_92))) (let ((.def_97 (not .def_96))) (let ((.def_98 (and .def_97 .def_90))) (let ((.def_99 (and A2 A15))) (let ((.def_100 (= .def_26 A15))) (let ((.def_101 (= .def_100 .def_99))) (let ((.def_102 (not .def_101))) (let ((.def_103 (or A25 A5))) (let ((.def_104 (not .def_103))) (let ((.def_105 (= A3 A24))) (let ((.def_106 (or .def_105 .def_104))) (let ((.def_107 (not .def_106))) (let ((.def_108 (and .def_107 .def_102))) (let ((.def_109 (or .def_108 .def_98))) (let ((.def_110 (not .def_109))) (let ((.def_111 (and .def_110 .def_83))) (let ((.def_112 (not .def_111))) (let ((.def_113 (or .def_112 .def_54))) .def_113)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
