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
(assert (let ((.def_0 (or A15 A9))) (let ((.def_1 (not A13))) (let ((.def_2 (not A2))) (let ((.def_3 (= .def_2 .def_1))) (let ((.def_4 (and .def_3 .def_0))) (let ((.def_5 (not A20))) (let ((.def_6 (not A22))) (let ((.def_7 (and .def_6 .def_5))) (let ((.def_8 (not A23))) (let ((.def_9 (and .def_1 .def_8))) (let ((.def_10 (or .def_9 .def_7))) (let ((.def_11 (or .def_10 .def_4))) (let ((.def_12 (not A11))) (let ((.def_13 (not A17))) (let ((.def_14 (and .def_13 .def_12))) (let ((.def_15 (and .def_13 A20))) (let ((.def_16 (or .def_15 .def_14))) (let ((.def_17 (not .def_16))) (let ((.def_18 (not A18))) (let ((.def_19 (and A1 .def_18))) (let ((.def_20 (not A21))) (let ((.def_21 (or .def_18 .def_20))) (let ((.def_22 (not .def_21))) (let ((.def_23 (or .def_22 .def_19))) (let ((.def_24 (and .def_23 .def_17))) (let ((.def_25 (not .def_24))) (let ((.def_26 (or .def_25 .def_11))) (let ((.def_27 (not A15))) (let ((.def_28 (= A12 .def_27))) (let ((.def_29 (not A10))) (let ((.def_30 (and A7 .def_29))) (let ((.def_31 (not .def_30))) (let ((.def_32 (and .def_31 .def_28))) (let ((.def_33 (not .def_32))) (let ((.def_34 (not A16))) (let ((.def_35 (and .def_34 A12))) (let ((.def_36 (not .def_35))) (let ((.def_37 (not A3))) (let ((.def_38 (not A6))) (let ((.def_39 (or .def_38 .def_37))) (let ((.def_40 (not .def_39))) (let ((.def_41 (and .def_40 .def_36))) (let ((.def_42 (not .def_41))) (let ((.def_43 (= .def_42 .def_33))) (let ((.def_44 (not .def_43))) (let ((.def_45 (= A4 A0))) (let ((.def_46 (not A5))) (let ((.def_47 (and .def_46 A8))) (let ((.def_48 (and .def_47 .def_45))) (let ((.def_49 (not .def_48))) (let ((.def_50 (not A19))) (let ((.def_51 (and A6 .def_50))) (let ((.def_52 (not .def_51))) (let ((.def_53 (or .def_46 .def_34))) (let ((.def_54 (not .def_53))) (let ((.def_55 (and .def_54 .def_52))) (let ((.def_56 (not .def_55))) (let ((.def_57 (or .def_56 .def_49))) (let ((.def_58 (not .def_57))) (let ((.def_59 (or .def_58 .def_44))) (let ((.def_60 (and .def_59 .def_26))) (let ((.def_61 (not .def_60))) (let ((.def_62 (or A12 A10))) (let ((.def_63 (not A24))) (let ((.def_64 (or A15 .def_63))) (let ((.def_65 (and .def_64 .def_62))) (let ((.def_66 (not .def_65))) (let ((.def_67 (or .def_27 A4))) (let ((.def_68 (not .def_67))) (let ((.def_69 (not A1))) (let ((.def_70 (and .def_38 .def_69))) (let ((.def_71 (not .def_70))) (let ((.def_72 (or .def_71 .def_68))) (let ((.def_73 (and .def_72 .def_66))) (let ((.def_74 (and A21 .def_1))) (let ((.def_75 (not .def_74))) (let ((.def_76 (not A12))) (let ((.def_77 (or .def_76 .def_1))) (let ((.def_78 (not .def_77))) (let ((.def_79 (and .def_78 .def_75))) (let ((.def_80 (not .def_79))) (let ((.def_81 (or A13 .def_34))) (let ((.def_82 (or .def_5 .def_63))) (let ((.def_83 (or .def_82 .def_81))) (let ((.def_84 (not .def_83))) (let ((.def_85 (or .def_84 .def_80))) (let ((.def_86 (not .def_85))) (let ((.def_87 (= .def_86 .def_73))) (let ((.def_88 (or .def_8 A13))) (let ((.def_89 (and .def_46 A4))) (let ((.def_90 (not .def_89))) (let ((.def_91 (and .def_90 .def_88))) (let ((.def_92 (and A13 A0))) (let ((.def_93 (not .def_92))) (let ((.def_94 (and .def_34 .def_1))) (let ((.def_95 (or .def_94 .def_93))) (let ((.def_96 (or .def_95 .def_91))) (let ((.def_97 (not A8))) (let ((.def_98 (or A6 .def_97))) (let ((.def_99 (not .def_98))) (let ((.def_100 (or .def_1 .def_69))) (let ((.def_101 (not .def_100))) (let ((.def_102 (and .def_101 .def_99))) (let ((.def_103 (or .def_20 .def_37))) (let ((.def_104 (and .def_50 .def_50))) (let ((.def_105 (= .def_104 .def_103))) (let ((.def_106 (not .def_105))) (let ((.def_107 (or .def_106 .def_102))) (let ((.def_108 (and .def_107 .def_96))) (let ((.def_109 (not .def_108))) (let ((.def_110 (and .def_109 .def_87))) (let ((.def_111 (and .def_110 .def_61))) (let ((.def_112 (not .def_111))) .def_112))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)