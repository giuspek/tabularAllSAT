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
(declare-fun A14 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (or A0 A23))) (let ((.def_1 (or A4 A16))) (let ((.def_2 (or .def_1 .def_0))) (let ((.def_3 (not A19))) (let ((.def_4 (and .def_3 A15))) (let ((.def_5 (not A24))) (let ((.def_6 (not A22))) (let ((.def_7 (and .def_6 .def_5))) (let ((.def_8 (not .def_7))) (let ((.def_9 (or .def_8 .def_4))) (let ((.def_10 (or .def_9 .def_2))) (let ((.def_11 (not .def_10))) (let ((.def_12 (or A2 A16))) (let ((.def_13 (not .def_12))) (let ((.def_14 (not A5))) (let ((.def_15 (and .def_14 A23))) (let ((.def_16 (or .def_15 .def_13))) (let ((.def_17 (not .def_16))) (let ((.def_18 (not A6))) (let ((.def_19 (= .def_18 A7))) (let ((.def_20 (and A13 A13))) (let ((.def_21 (not .def_20))) (let ((.def_22 (and .def_21 .def_19))) (let ((.def_23 (not .def_22))) (let ((.def_24 (and .def_23 .def_17))) (let ((.def_25 (and .def_24 .def_11))) (let ((.def_26 (not A16))) (let ((.def_27 (and .def_26 A19))) (let ((.def_28 (not A10))) (let ((.def_29 (or A16 .def_28))) (let ((.def_30 (not .def_29))) (let ((.def_31 (and .def_30 .def_27))) (let ((.def_32 (not .def_31))) (let ((.def_33 (not A7))) (let ((.def_34 (and .def_33 A4))) (let ((.def_35 (and A2 .def_6))) (let ((.def_36 (not .def_35))) (let ((.def_37 (and .def_36 .def_34))) (let ((.def_38 (not .def_37))) (let ((.def_39 (and .def_38 .def_32))) (let ((.def_40 (not .def_39))) (let ((.def_41 (or A4 .def_6))) (let ((.def_42 (and A16 A2))) (let ((.def_43 (not .def_42))) (let ((.def_44 (or .def_43 .def_41))) (let ((.def_45 (not A23))) (let ((.def_46 (or A8 .def_45))) (let ((.def_47 (and A2 .def_14))) (let ((.def_48 (or .def_47 .def_46))) (let ((.def_49 (or .def_48 .def_44))) (let ((.def_50 (not .def_49))) (let ((.def_51 (or .def_50 .def_40))) (let ((.def_52 (not .def_51))) (let ((.def_53 (and .def_52 .def_25))) (let ((.def_54 (not A8))) (let ((.def_55 (or A3 .def_54))) (let ((.def_56 (or .def_5 A14))) (let ((.def_57 (or .def_56 .def_55))) (let ((.def_58 (not .def_57))) (let ((.def_59 (not A18))) (let ((.def_60 (and A21 .def_59))) (let ((.def_61 (not .def_60))) (let ((.def_62 (or .def_6 A3))) (let ((.def_63 (= .def_62 .def_61))) (let ((.def_64 (or .def_63 .def_58))) (let ((.def_65 (not .def_64))) (let ((.def_66 (not A4))) (let ((.def_67 (not A1))) (let ((.def_68 (= .def_67 .def_66))) (let ((.def_69 (and A5 A13))) (let ((.def_70 (or .def_69 .def_68))) (let ((.def_71 (not .def_70))) (let ((.def_72 (and .def_5 A21))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and A1 A22))) (let ((.def_75 (and .def_74 .def_73))) (let ((.def_76 (not .def_75))) (let ((.def_77 (= .def_76 .def_71))) (let ((.def_78 (not .def_77))) (let ((.def_79 (or .def_78 .def_65))) (let ((.def_80 (not .def_79))) (let ((.def_81 (and A15 A18))) (let ((.def_82 (not .def_81))) (let ((.def_83 (or A20 .def_67))) (let ((.def_84 (not .def_83))) (let ((.def_85 (and .def_84 .def_82))) (let ((.def_86 (not A11))) (let ((.def_87 (or A21 .def_86))) (let ((.def_88 (not .def_87))) (let ((.def_89 (not A12))) (let ((.def_90 (or A18 .def_89))) (let ((.def_91 (not .def_90))) (let ((.def_92 (or .def_91 .def_88))) (let ((.def_93 (and .def_92 .def_85))) (let ((.def_94 (not .def_93))) (let ((.def_95 (and A19 A19))) (let ((.def_96 (or A15 A7))) (let ((.def_97 (or .def_96 .def_95))) (let ((.def_98 (not .def_97))) (let ((.def_99 (not A20))) (let ((.def_100 (or .def_99 .def_14))) (let ((.def_101 (not .def_100))) (let ((.def_102 (not A14))) (let ((.def_103 (or .def_102 .def_14))) (let ((.def_104 (or .def_103 .def_101))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and .def_105 .def_98))) (let ((.def_107 (not .def_106))) (let ((.def_108 (or .def_107 .def_94))) (let ((.def_109 (and .def_108 .def_80))) (let ((.def_110 (or .def_109 .def_53))) .def_110))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
