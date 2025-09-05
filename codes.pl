
=head   
####complicate version
    my $flag_s = 1;
    my $flag_e = 0;
    my $nodesc = keys %{$nhr};
    foreach my $n ( sort { $$nhr{$a}{I} <=> $$nhr{$b}{I} } keys %{$nhr} ){
        my $d = $$nhr{$n}{D};
        print "\tN:$n D:$d\n" if $debug;    
        #print "Extract:$gene node:$n S:$anno{$gene}{$n}{S} E:$anno{$gene}{$n}{E} \n" if $debug;
        if( defined  $anno{$gene}{$n} ){ 
            if($d * $strand[$which] * $anno{$gene}{$n}{D} >= 0){
                my $as = $anno{$gene}{$n}{S};
                my $ae = $anno{$gene}{$n}{E};
                my $nl = $anno{$gene}{$n}{L};
                my ($rs,$re) =  defined $$nhr{$n}{E} ? @{$$nhr{$n}{E}} : (0,0);
                my $clen = 0;   
                if($rs == 0 and $re ==0  and $as == 0 and $ae == 0){
                    $olen_n += $nl;;
                    $clen =  $nl;
                    print "\t\tall full cover, len is node lenght $nl\n" if $debug;

                }else{
                    $olen_n += $nl - $rs - $re; 
                    if($rs + 1 > $nl - $ae or $as + 1 > $nl -$re){
                        print "\t\toverlap is zero $nl $rs $re $as $ae\n" if $debug;
                    }else{
                        my $overlap = min($nl-$ae, $nl - $re) - max($rs, $as);
                        my $or = $overlap / ($nl - $rs - $re);
                        print "\t\tpartial overlap $nl $rs $re $as $ae $overlap $or\n" if $debug;
                        $clen =  $overlap;
                    }
                    if($flag_sus ){
                        print "\t\tparitial have overhang, save it to susp\n" if $debug;
                        $susp{$n} = $$nhr{$n};
                    }
                }
                $olen += $clen;
                print "\t$n have correct direc, overlap: $clen ends: $rs $re $as $ae \n" if $debug;
            }else{
                print "\tdirectino mismatch $d * $strand[$which] * $anno{$gene}{$n}{D} \n" if $debug;
            }
            $flag_sus = 0 if ($olen_n > 5) ;
        }else{
            if($flag_sus){
                print "\tsuspent end: $n...\n" if $debug;
                $susp{$n} = $$nhr{$n};
            }else{
                print "\tdon't exists in anno, not a suspend..\n" if $debug;
            }
        }
    }   
    print "OVER return nodeslen $olen_n  len:$olen totalnode:$nodesc \n" if $debug;
    return ($olen,$olen_n, $nodesc ,\%susp);
}                                                                                                                       
                                                                                                                        
