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
(assert (let ((.def_0 (not A10))) (let ((.def_1 (and .def_0 A5))) (let ((.def_2 (not A8))) (let ((.def_3 (not A1))) (let ((.def_4 (= .def_3 .def_2))) (let ((.def_5 (and .def_4 .def_1))) (let ((.def_6 (not .def_5))) (let ((.def_7 (or A3 .def_0))) (let ((.def_8 (not .def_7))) (let ((.def_9 (not A0))) (let ((.def_10 (not A22))) (let ((.def_11 (or .def_10 .def_9))) (let ((.def_12 (not .def_11))) (let ((.def_13 (and .def_12 .def_8))) (let ((.def_14 (or .def_13 .def_6))) (let ((.def_15 (not .def_14))) (let ((.def_16 (not A5))) (let ((.def_17 (and A14 .def_16))) (let ((.def_18 (not .def_17))) (let ((.def_19 (not A17))) (let ((.def_20 (not A7))) (let ((.def_21 (or .def_20 .def_19))) (let ((.def_22 (not .def_21))) (let ((.def_23 (= .def_22 .def_18))) (let ((.def_24 (= A15 .def_20))) (let ((.def_25 (not A20))) (let ((.def_26 (not A21))) (let ((.def_27 (and .def_26 .def_25))) (let ((.def_28 (or .def_27 .def_24))) (let ((.def_29 (or .def_28 .def_23))) (let ((.def_30 (and .def_29 .def_15))) (let ((.def_31 (and A22 .def_9))) (let ((.def_32 (and .def_16 A23))) (let ((.def_33 (not .def_32))) (let ((.def_34 (and .def_33 .def_31))) (let ((.def_35 (not A9))) (let ((.def_36 (or A2 .def_35))) (let ((.def_37 (and A1 A5))) (let ((.def_38 (not .def_37))) (let ((.def_39 (= .def_38 .def_36))) (let ((.def_40 (or .def_39 .def_34))) (let ((.def_41 (not .def_40))) (let ((.def_42 (and A18 A10))) (let ((.def_43 (or .def_26 A9))) (let ((.def_44 (not .def_43))) (let ((.def_45 (= .def_44 .def_42))) (let ((.def_46 (not .def_45))) (let ((.def_47 (not A16))) (let ((.def_48 (not A18))) (let ((.def_49 (and .def_48 .def_47))) (let ((.def_50 (not A13))) (let ((.def_51 (not A19))) (let ((.def_52 (= .def_51 .def_50))) (let ((.def_53 (and .def_52 .def_49))) (let ((.def_54 (or .def_53 .def_46))) (let ((.def_55 (and .def_54 .def_41))) (let ((.def_56 (and .def_55 .def_30))) (let ((.def_57 (not .def_56))) (let ((.def_58 (not A23))) (let ((.def_59 (or .def_58 A9))) (let ((.def_60 (not .def_59))) (let ((.def_61 (not A6))) (let ((.def_62 (or .def_61 A17))) (let ((.def_63 (not .def_62))) (let ((.def_64 (or .def_63 .def_60))) (let ((.def_65 (not .def_64))) (let ((.def_66 (= .def_10 A3))) (let ((.def_67 (not .def_66))) (let ((.def_68 (and A5 A0))) (let ((.def_69 (or .def_68 .def_67))) (let ((.def_70 (and .def_69 .def_65))) (let ((.def_71 (not .def_70))) (let ((.def_72 (not A3))) (let ((.def_73 (and .def_72 A0))) (let ((.def_74 (or A12 .def_51))) (let ((.def_75 (not .def_74))) (let ((.def_76 (and .def_75 .def_73))) (let ((.def_77 (not .def_76))) (let ((.def_78 (not A24))) (let ((.def_79 (or .def_78 .def_2))) (let ((.def_80 (not .def_79))) (let ((.def_81 (or A18 .def_20))) (let ((.def_82 (not .def_81))) (let ((.def_83 (and .def_82 .def_80))) (let ((.def_84 (not .def_83))) (let ((.def_85 (and .def_84 .def_77))) (let ((.def_86 (not .def_85))) (let ((.def_87 (or .def_86 .def_71))) (let ((.def_88 (or A11 .def_19))) (let ((.def_89 (not .def_88))) (let ((.def_90 (and .def_51 .def_61))) (let ((.def_91 (and .def_90 .def_89))) (let ((.def_92 (not A2))) (let ((.def_93 (and .def_92 A2))) (let ((.def_94 (not .def_93))) (let ((.def_95 (= .def_35 .def_19))) (let ((.def_96 (or .def_95 .def_94))) (let ((.def_97 (and .def_96 .def_91))) (let ((.def_98 (not .def_97))) (let ((.def_99 (not .def_36))) (let ((.def_100 (or A1 A12))) (let ((.def_101 (not .def_100))) (let ((.def_102 (and .def_101 .def_99))) (let ((.def_103 (not .def_102))) (let ((.def_104 (not A4))) (let ((.def_105 (and A7 .def_104))) (let ((.def_106 (or .def_50 A1))) (let ((.def_107 (not .def_106))) (let ((.def_108 (or .def_107 .def_105))) (let ((.def_109 (or .def_108 .def_103))) (let ((.def_110 (or .def_109 .def_98))) (let ((.def_111 (= .def_110 .def_87))) (let ((.def_112 (and .def_111 .def_57))) (let ((.def_113 (not .def_112))) .def_113)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
