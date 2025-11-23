'use client';

import React, { useState, useRef, useEffect } from 'react';
import { Info } from 'lucide-react';

interface InfoTooltipProps {
    content: string;
    title?: string;
    className?: string;
}

export default function InfoTooltip({ content, title, className = '' }: InfoTooltipProps) {
    const [isVisible, setIsVisible] = useState(false);
    const [position, setPosition] = useState<'top' | 'bottom' | 'left' | 'right'>('top');
    const [tooltipStyle, setTooltipStyle] = useState<React.CSSProperties>({ opacity: 0 });
    const [isPositioned, setIsPositioned] = useState(false);
    const iconRef = useRef<HTMLDivElement>(null);
    const tooltipRef = useRef<HTMLDivElement>(null);

    // Check if className contains text-white to adjust icon color
    const isWhiteText = className.includes('text-white') || className.includes('text-white/');
    const iconColor = isWhiteText 
        ? 'text-white/80 hover:text-white' 
        : 'text-gray-400 hover:text-gray-600 dark:hover:text-gray-300';

    // Calculate tooltip position based on available space
    useEffect(() => {
        if (!isVisible || !iconRef.current) return;

        // Use requestAnimationFrame to ensure DOM is updated
        const updatePosition = () => {
            if (!iconRef.current || !tooltipRef.current) return;

            const iconRect = iconRef.current.getBoundingClientRect();
            const tooltipRect = tooltipRef.current.getBoundingClientRect();
            const viewportWidth = window.innerWidth;
            const viewportHeight = window.innerHeight;

            const spaceAbove = iconRect.top;
            const spaceBelow = viewportHeight - iconRect.bottom;
            const spaceLeft = iconRect.left;
            const spaceRight = viewportWidth - iconRect.right;

            // Use actual tooltip dimensions if available, otherwise estimate
            const tooltipWidth = tooltipRect.width > 0 ? tooltipRect.width : 320;
            const tooltipHeight = tooltipRect.height > 0 ? tooltipRect.height : 150;

            let newPosition: 'top' | 'bottom' | 'left' | 'right' = 'top';

            // Determine best position based on available space
            if (spaceBelow >= tooltipHeight + 30) {
                newPosition = 'bottom';
            } else if (spaceAbove >= tooltipHeight + 30) {
                newPosition = 'top';
            } else if (spaceRight >= tooltipWidth + 30) {
                newPosition = 'right';
            } else if (spaceLeft >= tooltipWidth + 30) {
                newPosition = 'left';
            } else {
                // Default to bottom if there's more space there
                newPosition = spaceBelow > spaceAbove ? 'bottom' : 'top';
            }

            // Calculate position based on chosen side
            const iconCenterX = iconRect.left + iconRect.width / 2;
            const iconCenterY = iconRect.top + iconRect.height / 2;
            
            let top = 0;
            let left = 0;

            if (newPosition === 'bottom') {
                top = iconRect.bottom + 8;
                left = iconCenterX - tooltipWidth / 2;
            } else if (newPosition === 'top') {
                top = iconRect.top - tooltipHeight - 8;
                left = iconCenterX - tooltipWidth / 2;
            } else if (newPosition === 'right') {
                top = iconCenterY - tooltipHeight / 2;
                left = iconRect.right + 8;
            } else if (newPosition === 'left') {
                top = iconCenterY - tooltipHeight / 2;
                left = iconRect.left - tooltipWidth - 8;
            }

            // Ensure tooltip stays within viewport with padding
            const padding = 10;
            top = Math.max(padding, Math.min(top, viewportHeight - tooltipHeight - padding));
            left = Math.max(padding, Math.min(left, viewportWidth - tooltipWidth - padding));

            setPosition(newPosition);
            setTooltipStyle({
                position: 'fixed',
                top: `${top}px`,
                left: `${left}px`,
                zIndex: 9999,
                opacity: 1,
            });
            setIsPositioned(true);
        };

        // Wait for next frame for tooltip to render
        requestAnimationFrame(() => {
            requestAnimationFrame(updatePosition);
        });
    }, [isVisible]);

    // Handle window resize and scroll
    useEffect(() => {
        if (!isVisible) {
            setIsPositioned(false);
            return;
        }
        
        const handleReposition = () => {
            if (iconRef.current && tooltipRef.current) {
                const updatePosition = () => {
                    if (!iconRef.current || !tooltipRef.current) return;
                    
                    const iconRect = iconRef.current.getBoundingClientRect();
                    const tooltipRect = tooltipRef.current.getBoundingClientRect();
                    const viewportWidth = window.innerWidth;
                    const viewportHeight = window.innerHeight;
                    
                    const tooltipWidth = tooltipRect.width > 0 ? tooltipRect.width : 320;
                    const tooltipHeight = tooltipRect.height > 0 ? tooltipRect.height : 150;
                    
                    const spaceAbove = iconRect.top;
                    const spaceBelow = viewportHeight - iconRect.bottom;
                    
                    let newPosition: 'top' | 'bottom' | 'left' | 'right' = spaceBelow > spaceAbove ? 'bottom' : 'top';
                    const iconCenterX = iconRect.left + iconRect.width / 2;
                    
                    let top = newPosition === 'bottom' ? iconRect.bottom + 8 : iconRect.top - tooltipHeight - 8;
                    let left = iconCenterX - tooltipWidth / 2;
                    
                    const padding = 10;
                    top = Math.max(padding, Math.min(top, viewportHeight - tooltipHeight - padding));
                    left = Math.max(padding, Math.min(left, viewportWidth - tooltipWidth - padding));
                    
                    setPosition(newPosition);
                    setTooltipStyle(prev => ({
                        ...prev,
                        top: `${top}px`,
                        left: `${left}px`,
                    }));
                };
                requestAnimationFrame(updatePosition);
            }
        };
        
        window.addEventListener('resize', handleReposition);
        window.addEventListener('scroll', handleReposition, true);
        return () => {
            window.removeEventListener('resize', handleReposition);
            window.removeEventListener('scroll', handleReposition, true);
            setIsPositioned(false);
        };
    }, [isVisible]);

    const getArrowClass = () => {
        switch (position) {
            case 'bottom':
                return 'absolute -top-1 left-1/2 transform -translate-x-1/2';
            case 'top':
                return 'absolute -bottom-1 left-1/2 transform -translate-x-1/2';
            case 'right':
                return 'absolute -left-1 top-1/2 transform -translate-y-1/2';
            case 'left':
                return 'absolute -right-1 top-1/2 transform -translate-y-1/2';
            default:
                return 'absolute -bottom-1 left-1/2 transform -translate-x-1/2';
        }
    };

    const getArrowTransform = () => {
        switch (position) {
            case 'top':
                return 'rotate-45';
            case 'bottom':
                return 'rotate-45';
            case 'right':
                return '-rotate-45';
            case 'left':
                return 'rotate-135';
            default:
                return 'rotate-45';
        }
    };

    const getArrowBorder = () => {
        switch (position) {
            case 'bottom':
                return 'border-t border-l border-gray-700';
            case 'top':
                return 'border-b border-r border-gray-700';
            case 'right':
                return 'border-t border-l border-gray-700';
            case 'left':
                return 'border-b border-r border-gray-700';
            default:
                return 'border-b border-r border-gray-700';
        }
    };

    return (
        <>
            <div 
                ref={iconRef}
                className={`relative inline-flex items-center ${isWhiteText ? '' : className}`}
                onMouseEnter={() => setIsVisible(true)}
                onMouseLeave={() => setIsVisible(false)}
            >
                <Info className={`w-4 h-4 ${iconColor} cursor-help transition-colors flex-shrink-0`} />
            </div>
            
            {isVisible && (
                <div
                    ref={tooltipRef}
                    className="fixed w-64 sm:w-80 pointer-events-auto transition-opacity duration-150"
                    style={tooltipStyle}
                    onMouseEnter={() => setIsVisible(true)}
                    onMouseLeave={() => setIsVisible(false)}
                >
                    <div className="bg-gray-900 dark:bg-gray-800 text-white text-xs rounded-lg shadow-2xl p-3 border border-gray-700 relative">
                        {title && (
                            <div className="font-semibold mb-1 text-sm">{title}</div>
                        )}
                        <div className="text-gray-200 leading-relaxed whitespace-normal">
                            {content}
                        </div>
                        {/* Arrow */}
                        <div className={getArrowClass()}>
                            <div 
                                className={`w-2 h-2 bg-gray-900 dark:bg-gray-800 ${getArrowBorder()} transform ${getArrowTransform()}`}
                            />
                        </div>
                    </div>
                </div>
            )}
        </>
    );
}

