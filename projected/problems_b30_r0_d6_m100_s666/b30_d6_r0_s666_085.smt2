(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A12 () Bool)
(declare-fun A14 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(declare-fun A25 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(assert (let ((.def_0 (or A24 A10))) (let ((.def_1 (not .def_0))) (let ((.def_2 (not A17))) (let ((.def_3 (and A19 .def_2))) (let ((.def_4 (not .def_3))) (let ((.def_5 (and .def_4 .def_1))) (let ((.def_6 (not .def_5))) (let ((.def_7 (and A7 A25))) (let ((.def_8 (or A22 A24))) (let ((.def_9 (and .def_8 .def_7))) (let ((.def_10 (= .def_9 .def_6))) (let ((.def_11 (not A6))) (let ((.def_12 (not A20))) (let ((.def_13 (and .def_12 .def_11))) (let ((.def_14 (not A14))) (let ((.def_15 (not A1))) (let ((.def_16 (or .def_15 .def_14))) (let ((.def_17 (not .def_16))) (let ((.def_18 (= .def_17 .def_13))) (let ((.def_19 (or A9 A12))) (let ((.def_20 (not .def_19))) (let ((.def_21 (not A5))) (let ((.def_22 (and A19 .def_21))) (let ((.def_23 (or .def_22 .def_20))) (let ((.def_24 (and .def_23 .def_18))) (let ((.def_25 (not .def_24))) (let ((.def_26 (and .def_25 .def_10))) (let ((.def_27 (and A0 A16))) (let ((.def_28 (not .def_27))) (let ((.def_29 (not A3))) (let ((.def_30 (or .def_29 A3))) (let ((.def_31 (not .def_30))) (let ((.def_32 (and .def_31 .def_28))) (let ((.def_33 (not .def_32))) (let ((.def_34 (or .def_12 A3))) (let ((.def_35 (or A8 A9))) (let ((.def_36 (or .def_35 .def_34))) (let ((.def_37 (and .def_36 .def_33))) (let ((.def_38 (not .def_37))) (let ((.def_39 (or A21 A19))) (let ((.def_40 (not .def_39))) (let ((.def_41 (not A24))) (let ((.def_42 (or A7 .def_41))) (let ((.def_43 (not .def_42))) (let ((.def_44 (or .def_43 .def_40))) (let ((.def_45 (not A21))) (let ((.def_46 (or .def_2 .def_45))) (let ((.def_47 (or A25 A12))) (let ((.def_48 (not .def_47))) (let ((.def_49 (or .def_48 .def_46))) (let ((.def_50 (or .def_49 .def_44))) (let ((.def_51 (not .def_50))) (let ((.def_52 (or .def_51 .def_38))) (let ((.def_53 (and .def_52 .def_26))) (let ((.def_54 (not A27))) (let ((.def_55 (or A27 .def_54))) (let ((.def_56 (not .def_55))) (let ((.def_57 (and A24 A14))) (let ((.def_58 (not .def_57))) (let ((.def_59 (and .def_58 .def_56))) (let ((.def_60 (or A24 A23))) (let ((.def_61 (and A6 .def_45))) (let ((.def_62 (not .def_61))) (let ((.def_63 (or .def_62 .def_60))) (let ((.def_64 (and .def_63 .def_59))) (let ((.def_65 (not .def_64))) (let ((.def_66 (not A12))) (let ((.def_67 (or .def_66 A9))) (let ((.def_68 (not .def_67))) (let ((.def_69 (not A8))) (let ((.def_70 (and A20 .def_69))) (let ((.def_71 (not .def_70))) (let ((.def_72 (or .def_71 .def_68))) (let ((.def_73 (or A2 A17))) (let ((.def_74 (not .def_73))) (let ((.def_75 (not A9))) (let ((.def_76 (and .def_75 .def_21))) (let ((.def_77 (not .def_76))) (let ((.def_78 (or .def_77 .def_74))) (let ((.def_79 (not .def_78))) (let ((.def_80 (and .def_79 .def_72))) (let ((.def_81 (not .def_80))) (let ((.def_82 (and .def_81 .def_65))) (let ((.def_83 (not .def_82))) (let ((.def_84 (or .def_75 A6))) (let ((.def_85 (not .def_84))) (let ((.def_86 (and A23 .def_11))) (let ((.def_87 (or .def_86 .def_85))) (let ((.def_88 (not .def_87))) (let ((.def_89 (and .def_14 .def_14))) (let ((.def_90 (not .def_89))) (let ((.def_91 (or A25 A23))) (let ((.def_92 (not .def_91))) (let ((.def_93 (or .def_92 .def_90))) (let ((.def_94 (not .def_93))) (let ((.def_95 (and .def_94 .def_88))) (let ((.def_96 (not .def_95))) (let ((.def_97 (and A19 A14))) (let ((.def_98 (not A10))) (let ((.def_99 (or A6 .def_98))) (let ((.def_100 (not .def_99))) (let ((.def_101 (= .def_100 .def_97))) (let ((.def_102 (not .def_101))) (let ((.def_103 (not A28))) (let ((.def_104 (= .def_103 .def_12))) (let ((.def_105 (not .def_104))) (let ((.def_106 (or .def_41 A10))) (let ((.def_107 (and .def_106 .def_105))) (let ((.def_108 (or .def_107 .def_102))) (let ((.def_109 (not .def_108))) (let ((.def_110 (or .def_109 .def_96))) (let ((.def_111 (or .def_110 .def_83))) (let ((.def_112 (or .def_111 .def_53))) (let ((.def_113 (not .def_112))) .def_113)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)